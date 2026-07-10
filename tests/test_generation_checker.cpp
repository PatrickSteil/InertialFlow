#include "external/doctest.h"
#include "generation_checker.hpp"

TEST_CASE("GenerationChecker: freshly constructed checker has nothing marked") {
  GenerationChecker<uint16_t> gc(5);
  for (std::size_t i = 0; i < 5; ++i) {
    CHECK_FALSE(gc.isMarked(i));
  }
}

TEST_CASE("GenerationChecker: mark() sets isMarked() for that index only") {
  GenerationChecker<uint16_t> gc(5);
  gc.mark(2);
  CHECK(gc.isMarked(2));
  CHECK_FALSE(gc.isMarked(0));
  CHECK_FALSE(gc.isMarked(1));
  CHECK_FALSE(gc.isMarked(3));
  CHECK_FALSE(gc.isMarked(4));
}

TEST_CASE("GenerationChecker: reset() clears all previous marks") {
  GenerationChecker<uint16_t> gc(4);
  gc.mark(0);
  gc.mark(3);
  gc.reset();
  CHECK_FALSE(gc.isMarked(0));
  CHECK_FALSE(gc.isMarked(3));

  // New marks after reset should still work.
  gc.mark(1);
  CHECK(gc.isMarked(1));
}

TEST_CASE("GenerationChecker: resize() clears marks and shrinks/grows capacity") {
  GenerationChecker<uint16_t> gc(3);
  gc.mark(0);
  gc.mark(1);

  gc.resize(6);
  CHECK(gc.isValid(5));
  CHECK_FALSE(gc.isValid(6));
  // All marks should be cleared by resize.
  for (std::size_t i = 0; i < 6; ++i) {
    CHECK_FALSE(gc.isMarked(i));
  }
}

TEST_CASE("GenerationChecker: isValid() respects bounds") {
  GenerationChecker<uint16_t> gc(3);
  CHECK(gc.isValid(0));
  CHECK(gc.isValid(2));
  CHECK_FALSE(gc.isValid(3));
  CHECK_FALSE(gc.isValid(100));
}

TEST_CASE("GenerationChecker: generation counter overflow still behaves correctly") {
  // Use a tiny generation type so we can cheaply wrap it around.
  GenerationChecker<uint8_t> gc(4);

  gc.mark(0);
  CHECK(gc.isMarked(0));

  // Force enough reset() calls to wrap the 8-bit generation counter at
  // least once (255 resets moves generation from 1 up through 256 -> wraps
  // to 0 -> reset() detects this and bumps to 1, clearing the seen array).
  for (int i = 0; i < 300; ++i) {
    gc.reset();
  }

  // After wraparound, nothing should spuriously read as marked.
  for (std::size_t i = 0; i < 4; ++i) {
    CHECK_FALSE(gc.isMarked(i));
  }

  gc.mark(2);
  CHECK(gc.isMarked(2));
  CHECK_FALSE(gc.isMarked(0));
}
