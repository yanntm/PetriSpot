/*
 * Log.h
 *
 * Timestamped INFO logging to standard output, shared by every module.
 */
#ifndef PETRI_CORE_LOG_H_
#define PETRI_CORE_LOG_H_

#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>

namespace petri
{

/** Print "[YYYY-MM-DD HH:MM:SS] [INFO   ] message" on standard output. */
inline void writeToLog (const std::string &message)
{
  auto now = std::chrono::system_clock::now ();
  std::time_t now_c = std::chrono::system_clock::to_time_t (now);
  struct std::tm *parts = std::localtime (&now_c);

  std::cout << "[" << (parts->tm_year + 1900) << "-" << std::setw (2)
      << std::setfill ('0') << (parts->tm_mon + 1) << "-" << std::setw (2)
      << std::setfill ('0') << parts->tm_mday << " " << std::setw (2)
      << std::setfill ('0') << parts->tm_hour << ":" << std::setw (2)
      << std::setfill ('0') << parts->tm_min << ":" << std::setw (2)
      << std::setfill ('0') << parts->tm_sec << "] " << "[INFO   ] "
      << message << std::endl;
}

} // namespace petri

#endif /* PETRI_CORE_LOG_H_ */
