#include <SDL.h>
#include <algorithm>
#include <assert.h>
#include <chrono>
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
  glm::vec3 min, max;
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

struct BVHPrimitiveInfo {
  BVHPrimitiveInfo() {
  }
  BVHPrimitiveInfo(size_t primitiveNumber, const AABB &bounds)
      : primitiveNumber(primitiveNumber), bounds(bounds), centroid(.5f * bounds.min + .5f * bounds.max) {
  }
  size_t primitiveNumber;
  AABB bounds;
  glm::vec3 centroid;
};

struct Object {
  AABB aabb;
};

glm::vec3 aabbCenter(AABB bounds) {
  return 0.5f * bounds.min + 0.5f * bounds.max;
}

uint32_t MaximumExtent(AABB aabb) {
  glm::vec3 d = aabb.max - aabb.min;
  if (d.x > d.y && d.x > d.z)
    return 0;
  else if (d.y > d.z)
    return 1;
  else
    return 2;
}

struct Node {
  enum Type { LEAF, NODE };
  Node *left, *right;
  AABB aabb;
  Type type;
};

uint32_t partitionObjects(Object *objects, const uint32_t nrObjects) {
  AABB centerBounds = {{FLT_MAX, FLT_MAX, FLT_MAX}, {-FLT_MAX, -FLT_MAX, -FLT_MAX}};
  for (uint32_t i = 0; i < nrObjects; i++)
    centerBounds = Union(centerBounds, aabbCenter(objects[i].aabb));
  int currentAxis = MaximumExtent(centerBounds);

  int mid = nrObjects >> 1;
  float axisMiddle = (centerBounds.min[currentAxis] + centerBounds.max[currentAxis]) / 2;

  Object *middlePtr = std::partition(objects, objects + nrObjects, [currentAxis, axisMiddle](const Object &object) {
    return aabbCenter(object.aabb)[currentAxis] < axisMiddle;
  });

  mid = uint32_t(middlePtr - objects);
  if (mid != 0 && mid != nrObjects) {
    return mid;
  }
  std::sort(objects, objects + nrObjects, [currentAxis](const Object &a, const Object &b) {
    return aabbCenter(a.aabb)[currentAxis] < aabbCenter(b.aabb)[currentAxis];
  });
  return nrObjects >> 1;
}

#ifdef _WIN32
#include <intrin.h>
uint32_t clz(uint32_t value) {
  unsigned long leadingZeros = 0;
  unsigned char res = _BitScanReverse(&leadingZeros, value);
  assert(res);
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

uint32_t findSplit(uint32_t *sortedMortonCodes, uint32_t first, uint32_t last) {
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
#include <algorithm>
#include <thread>

void orderBlockIncreasing(uint32_t *elements, uint32_t begin, uint32_t length) {
  uint32_t half = length;
  uint32_t outerEnd = begin + length;
  while ((half >>= 1) > 0) {
    uint32_t outerStepSize = half << 1;
    for (uint32_t i = begin; i < outerEnd; i += outerStepSize) {
      const uint32_t innerEnd = i + half;
      for (uint32_t j = i; j < innerEnd; j++) {
        if (elements[j] > elements[j + half])
          std::swap(elements[j], elements[j + half]);
      }
    }
  }
}

void orderBlockDecreasing(uint32_t *elements, uint32_t begin, uint32_t length) {
  uint32_t half = length;
  const uint32_t outerEnd = begin + length;
  while ((half >>= 1) > 0) {
    const uint32_t outerStepSize = half << 1;
    for (uint32_t i = begin; i < outerEnd; i += outerStepSize) {
      const uint32_t innerEnd = i + half;
      for (uint32_t j = i; j < innerEnd; j++) {
        if (elements[j] < elements[j + half])
          std::swap(elements[j], elements[j + half]);
      }
    }
  }
}

void sort(uint32_t *elements, uint32_t start, uint32_t end, uint32_t blockSize) {
  for (uint32_t i = start; i < end; i += blockSize * 2) {
    orderBlockIncreasing(elements, i, blockSize);
    orderBlockDecreasing(elements, i + blockSize, blockSize);
  }
}

void bitonicSort(uint32_t *elements, uint32_t count) {
  const uint32_t threadCount = std::thread::hardware_concurrency();
  void *stackMemory = _alloca(sizeof(std::thread) * threadCount);
  std::thread threads[32];

  for (uint32_t blockSize = 2; blockSize <= count; blockSize *= 2) {
    uint32_t nrOfBlocks = count / (blockSize * 2);
    uint32_t blockPerThread = std::max(nrOfBlocks / threadCount, 1u);
    uint32_t currentThreadCount = nrOfBlocks / blockPerThread;

    for (uint32_t threadIndex = 0; threadIndex < currentThreadCount; threadIndex++) {
      uint32_t start = threadIndex * blockPerThread * blockSize * 2;
      uint32_t end = (threadIndex + 1) * blockPerThread * blockSize * 2;
      threads[threadIndex] = std::thread(sort, elements, start, end, blockSize);
    }
     for (uint32_t threadIndex = 0; threadIndex < currentThreadCount; threadIndex++) {
      threads[threadIndex].join();
    }
  }
  orderBlockIncreasing(elements, 0, count);
}

// void sort(int arr[], int n, int order) {
//  bitonicSort(arr, 0, n, order);
//}
Node *TopDownBVTreeMorton(uint32_t *sortedMortonCodes, uint32_t *sortedObjectsIds, uint32_t first, uint32_t last) {
  if (first == last)
    return nullptr;

  uint32_t split = findSplit(sortedMortonCodes, first, last);
  Node *childA = TopDownBVTreeMorton(sortedMortonCodes, sortedObjectsIds, first, split);
  Node *childB = TopDownBVTreeMorton(sortedMortonCodes, sortedObjectsIds, first, split);

  Node *resNode = new Node();
  resNode->left = childA;
  resNode->left = childB;
  return resNode;
}

void TopDownBVTree(Node *node, Object *objects, const uint32_t nrObjects) {
  assert(nrObjects);
  AABB aabb = {{FLT_MAX, FLT_MAX, FLT_MAX}, {-FLT_MAX, -FLT_MAX, -FLT_MAX}};
  for (uint32_t i = 0; i < nrObjects; i++) {
    aabb = Union(aabb, objects[i].aabb);
  }
  node->aabb = aabb;
  if (nrObjects == 1) {
    node->type = Node::LEAF;
  } else {
    node->type = Node::NODE;
    uint32_t k = partitionObjects(objects, nrObjects);
    node->left = new Node{};
    node->right = new Node{};

    TopDownBVTree(node->left, objects, k);
    TopDownBVTree(node->right, objects + k, nrObjects - k);
  }
}

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

void renderNode(SDL_Renderer *renderer, Node *node) {
  SDL_Rect rec = getRect(node->aabb);
  if (node->type == Node::LEAF) {
    SDL_SetRenderDrawColor(renderer, 0, 150, 0, 255);
    SDL_RenderFillRect(renderer, &rec);
    SDL_SetRenderDrawColor(renderer, 150, 0, 0, 255);
    SDL_RenderDrawRect(renderer, &rec);
  } else {
    SDL_SetRenderDrawColor(renderer, 150, 0, 0, 255);
    SDL_RenderDrawRect(renderer, &rec);

    renderNode(renderer, node->left);
    renderNode(renderer, node->right);
  }
}

void renderTree(Node *tree);

namespace timer {
typedef std::chrono::time_point<std::chrono::high_resolution_clock> Time;

inline Time now() {
  return std::chrono::high_resolution_clock::now();
}
inline uint64_t durationInMs(const Time &start, const Time &end) {
  return uint64_t(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
}
} // namespace timer

int main() {
  int count = 16777216;
  std::vector<uint32_t> e, ee, ek, ekk;
  for (int i = 0; i < count; i++)
    e.push_back(rand());
  //e = {3, 1, 5, 2, 6, 8, 9, 10};
  {
    uint64_t totalTime = 0;
    for (int i = 0; i < 10; i++) {
      ee = e;
      auto start = timer::now();
      std::sort(ee.begin(), ee.end());
      auto end = timer::now();
      totalTime = timer::durationInMs(start, end);
      ek = ee;
    }
    std::cout << totalTime << std::endl;
  }
  {
    uint64_t totalTime = 0;
    for (int i = 0; i < 10; i++) {
      ee = e;
      auto start = timer::now();
      bitonicSort(ee.data(), uint32_t(ee.size()));
      auto end = timer::now();
      totalTime = timer::durationInMs(start, end);
      ekk = ee;
    }
    std::cout << totalTime << std::endl;
  }

  for (int i = 0; i < count; i++)
    if (ek[i] != ekk[i])
      printf("here we are\n");

  // if (true)
  //  return 0;

  //// uint32_t t = ctz(8);
  // Node *tree = new Node();
  // std::vector<Object> objects;

  // objects.push_back({AABB{{1, 1, 1}, {5, 5, 5}}});
  // objects.push_back({AABB{{49, 49, 49}, {51, 51, 51}}});
  // objects.push_back({AABB{{70, 70, 70}, {75, 75, 75}}});
  // objects.push_back({AABB{{55, 33, 11}, {60, 40, 11}}});

  // TopDownBVTree(tree, objects.data(), uint32_t(objects.size()));

  // renderTree(tree);

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
    renderNode(renderer, tree);
    SDL_RenderPresent(renderer);
  }
  SDL_DestroyRenderer(renderer);
  SDL_DestroyWindow(window);
  SDL_Quit();
}