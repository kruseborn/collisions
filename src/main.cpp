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

bool TestAABBAABB(AABB a, AABB b) {
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

// struct Node {
//  enum Type { LEAF, NODE };
//  Node *left, *right;
//  AABB aabb;
//  Type type;
//};
//
// uint32_t partitionObjects(Object *objects, const uint32_t nrObjects) {
//  AABB centerBounds = {{FLT_MAX, FLT_MAX, FLT_MAX}, {-FLT_MAX, -FLT_MAX, -FLT_MAX}};
//  for (uint32_t i = 0; i < nrObjects; i++)
//    centerBounds = Union(centerBounds, aabbCenter(objects[i].aabb));
//  int currentAxis = MaximumExtent(centerBounds);
//
//  int mid = nrObjects >> 1;
//  float axisMiddle = (centerBounds.min[currentAxis] + centerBounds.max[currentAxis]) / 2;
//
//  Object *middlePtr = std::partition(objects, objects + nrObjects, [currentAxis, axisMiddle](const Object &object) {
//    return aabbCenter(object.aabb)[currentAxis] < axisMiddle;
//  });
//
//  mid = uint32_t(middlePtr - objects);
//  if (mid != 0 && mid != nrObjects) {
//    return mid;
//  }
//  std::sort(objects, objects + nrObjects, [currentAxis](const Object &a, const Object &b) {
//    return aabbCenter(a.aabb)[currentAxis] < aabbCenter(b.aabb)[currentAxis];
//  });
//  return nrObjects >> 1;
//}

#ifdef _WIN32
#include <intrin.h>

#define SCREEN_WIDTH 1000
#define SCREEN_HEIGHT 1000

SDL_Rect getRect(AABB aabb) {
  SDL_Rect rec = {int(aabb.min.x), int(aabb.min.y), int(aabb.max.x - aabb.min.x), int(aabb.max.y - aabb.min.y)};
  rec.x *= 10;
  rec.y *= 10;
  rec.w *= 10;
  rec.h *= 10;
  return rec;
}

// void renderNode(SDL_Renderer *renderer, Node *node) {
//  SDL_Rect rec = getRect(node->aabb);
//  if (node->type == Node::LEAF) {
//    SDL_SetRenderDrawColor(renderer, 0, 150, 0, 255);
//    SDL_RenderFillRect(renderer, &rec);
//    SDL_SetRenderDrawColor(renderer, 150, 0, 0, 255);
//    SDL_RenderDrawRect(renderer, &rec);
//  } else {
//    SDL_SetRenderDrawColor(renderer, 150, 0, 0, 255);
//    SDL_RenderDrawRect(renderer, &rec);
//
//    renderNode(renderer, node->left);
//    renderNode(renderer, node->right);
//  }
//}

// void renderTree(Node *tree);

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

glm::uvec2 determineRange(std::vector<uint32_t> &sortedMortonCodes, uint32_t count, uint32_t index) {
  int32_t commonPrefixLeft = commonPrefix(sortedMortonCodes, count, index, index - 1);
  int32_t commonPrefixRight = commonPrefix(sortedMortonCodes, count, index, index + 1);

  // find direction of the range
  int32_t d = commonPrefixLeft > commonPrefixRight ? -1 : 1;
  int32_t commonPrefixMin = std::min(commonPrefixLeft, commonPrefixRight);

  // move to outside range
  int32_t rangeMax = 2;
  while (commonPrefix(sortedMortonCodes, count, index, index + d * rangeMax) > commonPrefixMin) {
    rangeMax *= 2;
  }

  // binary search to find range index
  int32_t rangeLength = 0;
  while (rangeMax >= 1) {
    bool larger = commonPrefix(sortedMortonCodes, count, index, index + (rangeLength + rangeMax) * d) > commonPrefixMin;
    rangeLength += larger ? rangeMax : 0;
    rangeMax = rangeMax >> 1;
  }

  int i = index;
  d = (commonPrefix(sortedMortonCodes, count, i, i + 1) - commonPrefix(sortedMortonCodes, count, i, i - 1)) > 0 ? 1
                                                                                                                : -1;
  int commonPrefixMin2 = commonPrefix(sortedMortonCodes, count, i, i - d);
  int l_max = 2;
  while (commonPrefix(sortedMortonCodes, count, i, i + d * l_max) > commonPrefixMin2) {
    l_max *= 2;
  }

  int l = 0;
  int t = l_max;
  do {
    t = (t + 1) >> 1; // exponential decrease
    if (commonPrefix(sortedMortonCodes, count, i, i + d * (l + t)) > commonPrefixMin) {
      l += t;
    }
  } while (t > 1);

  if (l != rangeLength)
    printf("here we are\n");
  int32_t j = index + (rangeLength * d);
  glm::uvec2 range = d > 0 ? glm::uvec2(index, j) : glm::uvec2(j, index);
  return range;
}
struct Node {
  uint32_t objectId;
  Node *childA, *childB;
  Node *parent;
  AABB aabb;
};

void updateAABBFromBottom(Node *leafNode) {
  Node *stack[32] = {};
  uint32_t top = 0;
  stack[top++] = leafNode->parent;
  AABB aabb = leafNode->aabb;
  while (top > 0) {
    Node *node = stack[--top];
    node->aabb = Union(node->aabb, aabb);
    aabb = node->aabb;
    if (node->parent != nullptr) {
      assert(top < 32);
      stack[top++] = leafNode->parent;
    }
  }
}

Node *generateHirarchy(std::vector<uint32_t> &sortedMortonCodes, std::vector<AABB> &sortedObjectsIds,
                       uint32_t count) {
  std::vector<Node> leafNodes(count);
  std::vector<Node> internalNodes(count - 1);

  for (int i = 0; i < count; i++)
    leafNodes[i].aabb = sortedObjectsIds[i];

  for (uint32_t i = 0; i < count - 1; i++) {
    glm::vec2 range = determineRange(sortedMortonCodes, count, i);
    uint32_t split = findSplit(sortedMortonCodes, range.x, range.y);
    Node *childA;
    if (split == range.x)
      childA = &leafNodes[split];
    else
      childA = &internalNodes[split];

    Node *childB;
    if (split + 1 == range.y)
      childB = &leafNodes[split + 1];
    else
      childB = &internalNodes[split + 1];

    internalNodes[i].childA = childA;
    internalNodes[i].childB = childB;
    childA->parent = &internalNodes[i];
    childB->parent = &internalNodes[i];
  }

  // update bvh
  for (int i = 0; i < count; i++) {
    updateAABBFromBottom(&leafNodes[i]);
  }

  return &internalNodes[0];
}

// void traverseIterative(Node *root, AABB &queryAABB) {
//  // Allocate traversal stack from thread-local memory,
//  // and push NULL to indicate that there are no postponed nodes.
//  Node* stack[64];
//  Node **stackPtr = stack;
//  *stackPtr++ = nullptr; // push
//
//  // Traverse nodes starting from the root.
//  Node *node = root;
//  do {
//    // Check each child node for overlap.
//    Node *childL = node->childA;
//    Node *childR = node->childB;
//    bool overlapL = (checkOverlap(queryAABB, bvh.getAABB(childL)));
//    bool overlapR = (checkOverlap(queryAABB, bvh.getAABB(childR)));
//
//    // Query overlaps a leaf node => report collision.
//    if (overlapL && bvh.isLeaf(childL))
//      list.add(queryObjectIdx, bvh.getObjectIdx(childL));
//
//    if (overlapR && bvh.isLeaf(childR))
//      list.add(queryObjectIdx, bvh.getObjectIdx(childR));
//
//    // Query overlaps an internal node => traverse.
//    bool traverseL = (overlapL && !bvh.isLeaf(childL));
//    bool traverseR = (overlapR && !bvh.isLeaf(childR));
//
//    if (!traverseL && !traverseR)
//      node = *--stackPtr; // pop
//    else {
//      node = (traverseL) ? childL : childR;
//      if (traverseL && traverseR)
//        *stackPtr++ = childR; // push
//    }
//  } while (node != NULL);
//}







int main() {
  const int count = 4;
  std::vector<glm::vec3> positions(count);
  std::vector<AABB> aabbs(count);
  std::vector<uint32_t> mortonCodes(count);
  for (int i = 0; i < count; i++) {
    if (i == 0)
      positions[i] = {0.1, 0.1, 0.1};
    if (i == 1)
      positions[i] = {0.1, 0.1, 0.9};
    if (i == 2)
      positions[i] = {0.1, 0.9, 0.9};
    if (i == 3)
      positions[i] = {0.9, 0.9, 0.9};

    //positions[i].x = (rand() / float(RAND_MAX));
    //positions[i].y = (rand() / float(RAND_MAX));
    //positions[i].z = (rand() / float(RAND_MAX));
    aabbs[i] = {{positions[i] - glm::vec3{0.05}}, {positions[i] + glm::vec3{0.05}}};
    aabbs[i].position = positions[i];
    mortonCodes[i] = morton3D(positions[i].x, positions[i].y, positions[i].z);
    aabbs[i].morton = mortonCodes[i];
  }
  radixSort(mortonCodes.data(), aabbs.data(), mortonCodes.size());
  Node *node = generateHirarchy(mortonCodes, aabbs, count);

  //uint32_t mortonsortedArray[8] = {0b00001, 0b00010, 0b00100, 0b00101, 0b10011, 0b11000, 0b11001, 0b11110};
  //uint32_t idssortedArray[8] = {0, 1, 2, 3, 4, 5, 6, 7};
  //std::vector<uint32_t> a = {mortonsortedArray, mortonsortedArray + 8};
  //std::vector<uint32_t> b = {idssortedArray, idssortedArray + 8};

  

  return 0;
}

void renderTree(Node *tree) {
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
    // renderNode(renderer, tree);
    SDL_RenderPresent(renderer);
  }
  SDL_DestroyRenderer(renderer);
  SDL_DestroyWindow(window);
  SDL_Quit();
}