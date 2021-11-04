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
  AABB aabb[4];
  int count;
};

struct Qbvh {
  std::vector<Node4> nodes;
};



std::vector<AABB> traverseIterative4(Qbvh *bvh, AABB &queryAABB) {
  std::vector<AABB> res;

  // Allocate traversal stack from thread-local memory,
  // and push NULL to indicate that there are no postponed nodes.
  Node4 *stack[64];
  Node4 **stackPtr = stack;
  *stackPtr++ = nullptr; // push
  *stackPtr++ = &bvh->nodes[0]; // push

  // Traverse nodes starting from the root.
  bool overlap[4] = {};
  bool shouldTraverse[4] = {};

  while (true) {
    Node4 *node = *--stackPtr; // pop
    if (node == nullptr)
      break;

    for (int i = 0; i < node->count; i++) {
      overlap[i] = checkOverlap(queryAABB, node->aabb[i]);
    }

    for (int i = 0; i < node->count; i++) {
      if (overlap[i] && node->leafs[i])
        res.push_back(node->leafs[i]->aabb);
    }

    for (int i = 0; i < node->count; i++) {
      if (overlap[i] && !node->leafs[i]) {
        *stackPtr++ = &bvh->nodes[node->children[i]]; // push
      }
    }
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
    result.aabb[i] = linearChildAabbs[i];
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
    //if (i == 0)
    //  positions[i] = {0.1, 0.1, 0.1};
    //if (i == 1)
    //  positions[i] = {0.5, 0.1, 0.1};
    //if (i == 2)
    //  positions[i] = {0.8, 0.1, 0.1};
    //if (i == 3)
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
    //renderNode(renderer, tree, 0);
    //renderNode4(renderer, qbvh, 0);

    int x, y;
    Uint32 buttons;

    SDL_PumpEvents(); // make sure we have the latest mouse state.

    buttons = SDL_GetMouseState(&x, &y);

    glm::vec3 positions = {x/float(SCREEN_WIDTH), y/float(SCREEN_HEIGHT), 0.1};
    AABB aabb = {{positions - glm::vec3{0.05}}, {positions + glm::vec3{0.05}}};
    //auto res = traverseIterative(tree, aabb);
    uint64_t start = SDL_GetPerformanceCounter();
    auto res = traverseIterative4(qbvh, aabb);
    uint64_t end = SDL_GetPerformanceCounter();
    
    float secondsElapsed = (end - start) / (float)SDL_GetPerformanceFrequency();
    
    printf("%f\n", secondsElapsed);

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