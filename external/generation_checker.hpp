#pragma once

#include <cassert>
#include <concepts>
#include <cstdint>
#include <vector>

// Class that tracks whether a vertex has been seen using a generational
// counter. This avoids having to clear the entire seen array by incrementing
// the generation counter.
template <std::integral GenerationType = std::uint16_t>
class GenerationChecker {
public:
  // Constructor initializes the seen array with the given size and sets the
  // generation to 1.
  explicit GenerationChecker(std::size_t size = 0)
      : seen(size, 0), generation(1) {}

  // Resizes the seen array and resets all values to zero.
  void resize(std::size_t size) {
    seen.assign(size, static_cast<GenerationType>(0));
    /* std::fill(seen.begin(), seen.end(), static_cast<GenerationType>(0)); */
    generation = 1;
  }

  // Increments the generation counter, effectively resetting all marks.
  // If the counter overflows, resets the entire seen array.
  void reset() {
    ++generation;
    if (generation == 0) {
      std::fill(seen.begin(), seen.end(), 0);
      ++generation;
    }
  }

  // Checks whether the given index is within bounds.
  inline bool isValid(std::size_t i) const { return i < seen.size(); }

  // Checks if a vertex at index `i` is marked as seen.
  inline bool isMarked(std::size_t i) const {
    assert(isValid(i));
    return seen[i] == generation;
  }

  // Marks a vertex at index `i` as seen by setting its value to the current
  // generation.
  inline void mark(std::size_t i) {
    assert(isValid(i));
    seen[i] = generation;
  }

  std::vector<GenerationType> seen;
  GenerationType generation;
};
