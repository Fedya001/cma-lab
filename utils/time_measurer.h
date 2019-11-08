#pragma once

#include <chrono>

template<class TimeUnit>
class TimeMeasurer {
 private:
  using Clock = std::chrono::high_resolution_clock;

 public:
  TimeMeasurer();
  void UpdateTimeStamp();

  [[nodiscard]] int64_t GetElapsedTime() const;
  template<class Function>
  [[nodiscard]] int64_t MeasureFunctionTime(Function function);

 private:
  Clock::time_point last_timestamp_;
};

template<class TimeUnit>
TimeMeasurer<TimeUnit>::TimeMeasurer()
    : last_timestamp_(Clock::now()) {}

template<class TimeUnit>
void TimeMeasurer<TimeUnit>::UpdateTimeStamp() {
  last_timestamp_ = Clock::now();
}

template<class TimeUnit>
int64_t TimeMeasurer<TimeUnit>::GetElapsedTime() const {
  return static_cast<int64_t>(std::chrono::duration_cast<TimeUnit>(Clock::now() - last_timestamp_).count());
}

template<class TimeUnit>
template<class Function>
int64_t TimeMeasurer<TimeUnit>::MeasureFunctionTime(Function function) {
  UpdateTimeStamp();
  function();
  return GetElapsedTime();
}
