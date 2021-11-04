#include <SDL.h>
#include <algorithm>
#include <assert.h>
#include <chrono>
#include <cmath>
#include <execution>
#include <float.h>
#include <glm/glm.hpp>
#include <inttypes.h>
#include <iostream>
#include <memory>
#include <stdio.h>
#include <string>
#include <vector>
#undef main
struct AABB {
  glm::vec3 min = {FLT_MAX, FLT_MAX, FLT_MAX};
  glm::vec3 max = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
  glm::vec3 position;
  uint32_t morton;
};

AABB Union(const AABB &a, const AABB &b) {
  AABB ret;
  ret.min = glm::min(a.min, b.min);
  ret.max = glm::max(a.max, b.max);
  return ret;
}

AABB Union(const AABB &a, glm::vec3 point) {
  AABB ret;
  ret.min = glm::min(a.min, point);
  ret.max = glm::max(a.max, point);
  return ret;
}

bool checkOverlap(AABB a, AABB b) {
  // Exit with no intersection if separated along an axis
  if (a.max[0] < b.min[0] || a.min[0] > b.max[0])
    return false;
  if (a.max[1] < b.min[1] || a.min[1] > b.max[1])
    return false;
  if (a.max[2] < b.min[2] || a.min[2] > b.max[2])
    return false;
  // Overlapping on all axes means AABBs are intersecting
  return true;
}

struct Object {
  AABB aabb;
};

glm::vec3 aabbCenter(AABB bounds) {
  return 0.5f * bounds.min + 0.5f * bounds.max;
}

#ifdef _WIN32
#include <intrin.h>

#define SCREEN_WIDTH 1000
#define SCREEN_HEIGHT 1000

namespace timer {
typedef std::chrono::time_point<std::chrono::high_resolution_clock> Time;

inline Time now() {
  return std::chrono::high_resolution_clock::now();
}
inline uint64_t durationInMs(const Time &start, const Time &end) {
  return uint64_t(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
}
} // namespace timer

static void radixSort(uint32_t *arr, AABB *aabbs, uint32_t size) {
  uint32_t *output = (uint32_t *)_malloca(sizeof(uint32_t) * size);
  AABB *output_aabbs = (AABB *)_malloca(sizeof(AABB) * size);

  uint32_t count[256];

  for (uint32_t s = 0; s < 32; s += 8) {
    memset((void *)count, 0, sizeof(count));
    for (uint32_t i = 0; i < size; i++) {
      ++count[(arr[i] >> s) & 0xff];
    }
    for (uint32_t i = 1; i < 256; i++) {
      count[i] += count[i - 1];
    }
    for (int32_t i = size - 1; i >= 0; i--) {
      uint32_t idx = (arr[i] >> s) & 0xff;
      output[--count[idx]] = arr[i];
      output_aabbs[count[idx]] = aabbs[i];
    }
    std::swap(arr, output);
    std::swap(aabbs, output_aabbs);
  }
  _freea(output);
  _freea(output_aabbs);
}

uint32_t clz(uint32_t value) {
  unsigned long leadingZeros = 0;
  leadingZeros = __lzcnt(value);
  return (leadingZeros);
}
#else
uint32_t ctz(uint32_t value) {
  return __builtin_clz(value);
}
#endif

// Expands a 10-bit integer into 30 bits
// by inserting 2 zeros after each bit.
static uint32_t expandBits(unsigned int v) {
  v = (v * 0x00010001u) & 0xFF0000FFu;
  v = (v * 0x00000101u) & 0x0F00F00Fu;
  v = (v * 0x00000011u) & 0xC30C30C3u;
  v = (v * 0x00000005u) & 0x49249249u;
  return v;
}

// Calculates a 30-bit Morton code for the
// given 3D point located within the unit cube [0,1].
static uint32_t morton3D(float x, float y, float z) {
  x = std::min(std::max(x * 1024.0f, 0.0f), 1023.0f);
  y = std::min(std::max(y * 1024.0f, 0.0f), 1023.0f);
  z = std::min(std::max(z * 1024.0f, 0.0f), 1023.0f);
  unsigned int xx = expandBits((unsigned int)x);
  unsigned int yy = expandBits((unsigned int)y);
  unsigned int zz = expandBits((unsigned int)z);
  return xx * 4 + yy * 2 + zz;
}

uint32_t findSplit(std::vector<uint32_t> &sortedMortonCodes, uint32_t first, uint32_t last) {
  uint32_t firstCode = sortedMortonCodes[first];
  uint32_t lastCode = sortedMortonCodes[last];
  if (first == lastCode)
    return (first + lastCode) >> 1;

  int commonPrefix = clz(firstCode ^ lastCode);
  int split = first; // initial guess
  int step = last - first;
  do {
    step = (step + 1) >> 1;           // exponential decrease
    uint32_t newSplit = split + step; // proposed new position

    if (newSplit < last) {
      unsigned int splitCode = sortedMortonCodes[newSplit];
      int splitPrefix = clz(firstCode ^ splitCode);
      if (splitPrefix > commonPrefix)
        split = newSplit;
    }
  } while (step > 1);
  return split;
}

inline int32_t commonPrefix(std::vector<uint32_t> &sortedMortonCodes, uint32_t count, uint32_t i, uint32_t j) {
  if (j < 0 || j >= count)
    return -1;
  assert(i != j);
  return sortedMortonCodes[i] ^ sortedMortonCodes[j];
}

struct Node {
  Node *childA, *childB;
  AABB aabb;
};

bool isLeaf(Node *node) {
  return node->childA == nullptr && node->childB == nullptr;
}

struct Ray {
  glm::vec3 O;
  glm::vec3 rD;
};

struct AABB4 {
  union {
    struct {
      float minX[4], minY[4], minZ[4];
      float maxX[4], maxY[4], maxZ[4];
    };
    struct {
      __m128 v_minX, v_minY, v_minZ, v_maxX, v_maxY, v_maxZ;
    };
  };
  AABB getAABB(uint32_t i) {
    return {{minX[i], minY[i], minZ[i]}, {maxX[i], maxY[i], maxZ[i]}};
  }
};

void intersection(Ray r, AABB4 b, bool res[4]) {
  __m128 rx = _mm_set1_ps(r.O.x);
  __m128 ry = _mm_set1_ps(r.O.y);
  __m128 rz = _mm_set1_ps(r.O.z);

  __m128 rdx = _mm_set1_ps(r.rD.x);
  __m128 rdy = _mm_set1_ps(r.rD.y);
  __m128 rdz = _mm_set1_ps(r.rD.z);

  __m128 zero = _mm_set1_ps(0.0f);

  // x plane
  __m128 t1 = _mm_mul_ps(_mm_sub_ps(b.v_minX, rx), rdx);
  __m128 t2 = _mm_mul_ps(_mm_sub_ps(b.v_maxX, rx), rdx);

  __m128 tmin = _mm_min_ps(t1, t2);
  __m128 tmax = _mm_max_ps(t1, t2);

  // y plane
  t1 = _mm_mul_ps(_mm_sub_ps(b.v_minY, ry), rdy);
  t2 = _mm_mul_ps(_mm_sub_ps(b.v_maxY, ry), rdy);

  tmin = _mm_max_ps(tmin, _mm_min_ps(t1, t2));
  tmax = _mm_min_ps(tmax, _mm_max_ps(t1, t2));

  // z plane
  t1 = _mm_mul_ps(_mm_sub_ps(b.v_minZ, rz), rdz);
  t2 = _mm_mul_ps(_mm_sub_ps(b.v_maxZ, rz), rdz);

  tmin = _mm_max_ps(tmin, _mm_min_ps(t1, t2));
  tmax = _mm_min_ps(tmax, _mm_max_ps(t1, t2));

  __m128 mask = _mm_cmple_ps(_mm_cmple_ps(tmax, tmin), zero);
  uint32_t *mask_ints = (uint32_t *)&mask;
  res[0] = mask_ints[0] > 0;
  res[1] = mask_ints[1] > 0;
  res[2] = mask_ints[2] > 0;
  res[3] = mask_ints[3] > 0;
}

bool intersection(Ray r, AABB b) {
  float tx1 = (b.min.x - r.O.x) * r.rD.x;
  float tx2 = (b.max.x - r.O.x) * r.rD.x;
  float tmin = std::min(tx1, tx2);
  float tmax = std::max(tx1, tx2);
  float ty1 = (b.min.y - r.O.y) * r.rD.y;
  float ty2 = (b.max.y - r.O.y) * r.rD.y;
  tmin = std::max(tmin, std::min(ty1, ty2));
  tmax = std::min(tmax, std::max(ty1, ty2));
  float tz1 = (b.min.z - r.O.z) * r.rD.z;
  float tz2 = (b.max.z - r.O.z) * r.rD.z;
  tmin = std::max(tmin, std::min(tz1, tz2));
  tmax = std::min(tmax, std::max(tz1, tz2));
  return tmax >= tmin && tmax >= 0;
}

// bool intersectionSIMD(Ray r, AABB b) {
//  __m128 t1 = _mm_mul_ps(_mm_sub_ps(node->bmin4, O4), rD4);
//  __m128 t2 = _mm_mul_ps(_mm_sub_ps(node->bmax4, O4), rD4);
//  __m128 vmax4 = _mm_max_ps(t1, t2), vmin4 = _mm_min_ps(t1, t2);
//  float *vmax = (float *)&vmax4, *vmin = (float *)&vmin4;
//  float tmax = min(vmax[0], min(vmax[1], vmax[2]));
//  float tmin = max(vmin[0], max(vmin[1], vmin[2]));
//  return tmax >= tmin && tmax >= 0;
//}

std::vector<AABB> traverseIterative(Node *root, AABB &queryAABB) {
  std::vector<AABB> res;

  // Allocate traversal stack from thread-local memory,
  // and push NULL to indicate that there are no postponed nodes.
  Node *stack[64];
  Node **stackPtr = stack;
  *stackPtr++ = NULL; // push

  // Traverse nodes starting from the root.
  Node *node = root;
  do {
    // Check each child node for overlap.
    Node *childL = node->childA;
    Node *childR = node->childB;
    bool overlapL = checkOverlap(queryAABB, childL->aabb);
    bool overlapR = checkOverlap(queryAABB, childR->aabb);

    // Query overlaps a leaf node => report collision.
    if (overlapL && isLeaf(childL))
      res.push_back(childL->aabb);

    if (overlapR && isLeaf(childR))
      res.push_back(childR->aabb);

    // Query overlaps an internal node => traverse.
    bool traverseL = (overlapL && !isLeaf(childL));
    bool traverseR = (overlapR && !isLeaf(childR));

    if (!traverseL && !traverseR)
      node = *--stackPtr; // pop
    else {
      node = (traverseL) ? childL : childR;
      if (traverseL && traverseR)
        *stackPtr++ = childR; // push
    }
  } while (node != NULL);

  return res;
}
struct Node4 {
  int32_t children[4];
  Node *leafs[4];
  AABB4 aabbs;
  int count;
};

struct Qbvh {
  std::vector<Node4> nodes;
};

// std::vector<AABB> traverseIterative4(Qbvh *bvh, AABB &queryAABB) {
//  std::vector<AABB> res;
//
//  Node4 *stack[64];
//  Node4 **stackPtr = stack;
//  *stackPtr++ = nullptr; // push
//
//  // Traverse nodes starting from the root.
//  // bool overlap[4] = {};
//  // bool shouldTraverse[4] = {};
//
//  Node4 *node = &bvh->nodes[0];
//  while (node != nullptr) {
//    for (int i = 0; i < node->count; i++) {
//      if (checkOverlap(queryAABB, node->aabb[i])) {
//        if (node->leafs[i])
//          res.push_back(node->leafs[i]->aabb);
//        else
//          *stackPtr++ = &bvh->nodes[node->children[i]]; // push
//      }
//    }
//    node = *--stackPtr; // pop
//  }
//  return res;
//}

std::vector<AABB> traverseIterative4(Qbvh *bvh, Ray ray) {
  std::vector<AABB> res;

  Node4 *stack[64];
  Node4 **stackPtr = stack;
  *stackPtr++ = nullptr; // push

  // Traverse nodes starting from the root.
  // bool overlap[4] = {};
  // bool shouldTraverse[4] = {};

  Node4 *node = &bvh->nodes[0];
  while (node != nullptr) {
    bool hit[4];
    intersection(ray, node->aabbs, hit);

    for (int i = 0; i < node->count; i++) {
      if (hit[i]) {
        if (node->leafs[i])
          res.push_back(node->leafs[i]->aabb);
        else
          *stackPtr++ = &bvh->nodes[node->children[i]]; // push
      }
    }
    node = *--stackPtr; // pop
  }
  return res;
}

std::vector<AABB> traverseIterative44(Qbvh *bvh, Ray ray) {
  std::vector<AABB> res;

  Node4 *stack[64];
  Node4 **stackPtr = stack;
  *stackPtr++ = nullptr; // push

  // Traverse nodes starting from the root.
  // bool overlap[4] = {};
  // bool shouldTraverse[4] = {};

  Node4 *node = &bvh->nodes[0];
  while (node != nullptr) {

    for (int i = 0; i < node->count; i++) {
      bool hit = intersection(ray, node->aabbs.getAABB(i));
      if (hit) {
        if (node->leafs[i])
          res.push_back(node->leafs[i]->aabb);
        else
          *stackPtr++ = &bvh->nodes[node->children[i]]; // push
      }
    }
    node = *--stackPtr; // pop
  }
  return res;
}

Node *generateHirarchy(std::vector<uint32_t> &sortedMortonCodes, std::vector<AABB> &sortedObjectAABBs, int32_t first,
                       int32_t last) {

  if (first == last) {
    Node *node = new Node{};
    node->aabb = sortedObjectAABBs[first];
    return node;
  }
  int split = findSplit(sortedMortonCodes, first, last);
  Node *childA = generateHirarchy(sortedMortonCodes, sortedObjectAABBs, first, split);
  Node *childB = generateHirarchy(sortedMortonCodes, sortedObjectAABBs, split + 1, last);

  Node *node = new Node{};
  node->aabb = Union(childA->aabb, childB->aabb);
  node->childA = childA;
  node->childB = childB;
  return node;
}

int constructQbvh(Qbvh *bvh, Node *node) {
  Node *linearChildIndices[4];
  AABB linearChildAabbs[4];
  int childCount = 0;

  {
    Node *currentNode = node->childA;
    if (!isLeaf(currentNode)) {
      linearChildIndices[childCount] = currentNode->childA;
      linearChildAabbs[childCount] = currentNode->childA->aabb;
      childCount++;
      linearChildIndices[childCount] = currentNode->childB;
      linearChildAabbs[childCount] = currentNode->childB->aabb;
      childCount++;
    } else {
      linearChildIndices[childCount] = currentNode;
      linearChildAabbs[childCount] = currentNode->aabb;
      childCount++;
    }
  }
  {
    Node *currentNode = node->childB;
    if (!isLeaf(currentNode)) {
      linearChildIndices[childCount] = currentNode->childA;
      linearChildAabbs[childCount] = currentNode->childA->aabb;
      childCount++;
      linearChildIndices[childCount] = currentNode->childB;
      linearChildAabbs[childCount] = currentNode->childB->aabb;
      childCount++;
    } else {
      linearChildIndices[childCount] = currentNode;
      linearChildAabbs[childCount] = currentNode->aabb;
      childCount++;
    }
  }
  Node4 result = {};
  result.count = childCount;
  for (int i = 0; i < childCount; i++) {
    result.aabbs.minX[i] = linearChildAabbs[i].min[0];
    result.aabbs.maxX[i] = linearChildAabbs[i].max[0];

    result.aabbs.minY[i] = linearChildAabbs[i].min[1];
    result.aabbs.maxY[i] = linearChildAabbs[i].max[1];

    result.aabbs.minZ[i] = linearChildAabbs[i].min[2];
    result.aabbs.maxZ[i] = linearChildAabbs[i].max[2];
  }
  size_t currentNodeIndex = bvh->nodes.size();
  bvh->nodes.push_back(result);
  int32_t resultIndex = static_cast<int32_t>(currentNodeIndex + 1);

  for (size_t i = 0; i < childCount; i++) {
    if (isLeaf(linearChildIndices[i])) { // is leaf
      bvh->nodes[currentNodeIndex].children[i] = -1;
      bvh->nodes[currentNodeIndex].leafs[i] = linearChildIndices[i];
    } else {
      bvh->nodes[currentNodeIndex].children[i] = resultIndex;
      resultIndex = constructQbvh(bvh, linearChildIndices[i]);
    }
  }
  return resultIndex;
}

void renderTree(Node *tree, Qbvh *qbvh);

int main() {
  const int count = 100000;
  std::vector<glm::vec3> positions(count);
  std::vector<AABB> aabbs(count);
  std::vector<uint32_t> mortonCodes(count);
  for (int i = 0; i < count; i++) {
    // if (i == 0)
    //  positions[i] = {0.5, 0.5, 0.1};
    // if (i == 1)
    //  positions[i] = {0.5, 0.1, 0.1};
    // if (i == 2)
    //  positions[i] = {0.8, 0.1, 0.1};
    // if (i == 3)
    //  positions[i] = {0.9, 0.1, 0.1};

    positions[i] = {rand() / float(RAND_MAX), rand() / float(RAND_MAX), 0.1};

    aabbs[i] = {{positions[i] - glm::vec3{0.01}}, {positions[i] + glm::vec3{0.01}}};
    aabbs[i].position = positions[i];
    mortonCodes[i] = morton3D(positions[i].x, positions[i].y, positions[i].z);
    aabbs[i].morton = mortonCodes[i];
  }
  radixSort(mortonCodes.data(), aabbs.data(), mortonCodes.size());
  Node *root = generateHirarchy(mortonCodes, aabbs, 0, count - 1);

  Qbvh qbvh;
  constructQbvh(&qbvh, root);

  renderTree(root, &qbvh);

  return 0;
}

SDL_Rect getRect(AABB aabb, float offsetY = 0) {
  aabb.min *= SCREEN_WIDTH;
  aabb.max *= SCREEN_WIDTH;

  aabb.min.y += offsetY;
  aabb.max.y += offsetY;
  SDL_Rect rec = {int(aabb.min.x), int(aabb.min.y), int(aabb.max.x - aabb.min.x), int(aabb.max.y - aabb.min.y)};
  return rec;
}

void renderAABB(SDL_Renderer *renderer, AABB aabb) {
  SDL_Rect rec = getRect(aabb);
  SDL_SetRenderDrawColor(renderer, 150, 0, 0, 50);
  SDL_RenderDrawRect(renderer, &rec);
}

void renderAABB_RED(SDL_Renderer *renderer, AABB aabb) {
  SDL_Rect rec = getRect(aabb);
  SDL_SetRenderDrawColor(renderer, 150, 150, 150, 255);
  SDL_RenderFillRect(renderer, &rec);
  SDL_SetRenderDrawColor(renderer, 0, 150, 0, 255);
  SDL_RenderDrawRect(renderer, &rec);
}
void renderNode4(SDL_Renderer *renderer, Qbvh *qbvh, uint32_t index) {
  Node4 node = qbvh->nodes[index];
  for (int i = 0; i < 4; i++) {
    if (node.children[i] > 0) {
      renderNode4(renderer, qbvh, node.children[i]);
    } else if (node.leafs[i]) {
      SDL_Rect rec = getRect(node.leafs[i]->aabb, 200.0f);
      SDL_SetRenderDrawColor(renderer, 0, 150, 150, 255);
      SDL_RenderFillRect(renderer, &rec);
      SDL_SetRenderDrawColor(renderer, 0, 0, 150, 255);
      SDL_RenderDrawRect(renderer, &rec);
    }
  }
}

void renderNode(SDL_Renderer *renderer, Node *node, int depth) {

  SDL_Rect rec = getRect(node->aabb);
  if (node->childA == nullptr) {
    SDL_SetRenderDrawColor(renderer, 0, 150, 150, 255);
    SDL_RenderFillRect(renderer, &rec);
    SDL_SetRenderDrawColor(renderer, 0, 0, 150, 255);
    SDL_RenderDrawRect(renderer, &rec);
  } else {
    renderNode(renderer, node->childA, depth + 1);
    renderNode(renderer, node->childB, depth + 1);
  }
}

void renderTree(Node *tree, Qbvh *qbvh) {
  if (SDL_Init(SDL_INIT_EVERYTHING) != 0) {
    exit(1);
  }
  SDL_Window *window = SDL_CreateWindow("hello_sdl2", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH,
                                        SCREEN_HEIGHT, SDL_WINDOW_SHOWN);
  if (window == NULL) {
    fprintf(stderr, "could not create window: %s\n", SDL_GetError());
    exit(1);
  }
  SDL_Surface *screenSurface = SDL_GetWindowSurface(window);
  SDL_Renderer *renderer = SDL_CreateRenderer(window, -1, 0);

  bool quit = false;
  while (!quit) {
    SDL_Event event;
    SDL_PollEvent(&event);
    switch (event.type) {
    case SDL_QUIT:
      quit = true;
      break;
    }

    SDL_SetRenderDrawColor(renderer, 242, 242, 242, 255);
    SDL_RenderClear(renderer);
    renderNode(renderer, tree, 0);
    // renderNode4(renderer, qbvh, 0);

    int x, y;
    Uint32 buttons;

    SDL_PumpEvents(); // make sure we have the latest mouse state.

    buttons = SDL_GetMouseState(&x, &y);

    glm::vec3 positions = {x / float(SCREEN_WIDTH), y / float(SCREEN_HEIGHT), 0.1};
    AABB aabb = {{positions - glm::vec3{0.05}}, {positions + glm::vec3{0.05}}};

    Ray ray;
    ray.O = positions;
    ray.O.z = 1;
    ray.rD = {1 / 1e-20, 1 / 1e-20, 1 / -1.0f};
    uint64_t start1 = SDL_GetPerformanceCounter();
    auto res = traverseIterative4(qbvh, ray);
    uint64_t end1 = SDL_GetPerformanceCounter();

    uint64_t start2 = SDL_GetPerformanceCounter();
    auto res2 = traverseIterative44(qbvh, ray);
    uint64_t end2 = SDL_GetPerformanceCounter();

    float secondsElapsed1 = (end1 - start1) / (float)SDL_GetPerformanceFrequency();
    float secondsElapsed2 = (end2 - start2) / (float)SDL_GetPerformanceFrequency();

    printf("%u %f %f %u %u\n", secondsElapsed1 < secondsElapsed2, secondsElapsed1, secondsElapsed2,
           (uint32_t)res.size(), (uint32_t)res2.size());

    for (int i = 0; i < res.size(); i++) {
      renderAABB_RED(renderer, res[i]);
    }
    renderAABB(renderer, aabb);

    SDL_RenderPresent(renderer);
  }
  SDL_DestroyRenderer(renderer);
  SDL_DestroyWindow(window);
  SDL_Quit();
}