// Taken and adapted from Ben Strasser

#pragma once

#include <chrono>
#include <iostream>
#include <string>

class StatusLog {
  StatusLog() = delete;
  StatusLog(const StatusLog &) = delete;
  StatusLog &operator=(const StatusLog &) = delete;

public:
  explicit StatusLog(const std::string &msg)
      : start_(std::chrono::steady_clock::now()) {
    std::cout << msg << " ... " << std::flush;
  }

  explicit StatusLog(const char *msg)
      : start_(std::chrono::steady_clock::now()) {
    std::cout << msg << " ... " << std::flush;
  }

  ~StatusLog() {
    const auto end = std::chrono::steady_clock::now();
    const auto ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start_)
            .count();

    std::cout << "done [" << ms << "ms]" << std::endl;
  }

private:
  std::chrono::steady_clock::time_point start_;
};
