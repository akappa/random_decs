//
// Created by andrea on 11/12/16.
//

#ifndef BROTLI_LOGGERS_H
#define BROTLI_LOGGERS_H

#include <ctime>
#include <cstdint>
#include <iostream>
#include <sstream>

enum class LoggerEvent {
    INFO, WARNING, ERROR
};

template <typename os>
os &operator<<(os &stream, LoggerEvent evt)
{
    switch (evt) {
        case LoggerEvent::WARNING: stream << "WARNING"; break;
        case LoggerEvent::ERROR: stream << "ERROR"; break;
        case LoggerEvent::INFO: stream << "INFO"; break;
        default: stream << "UNKNOWN";
    }
    return stream;
}

class Logger {
public:
    virtual void operator()(LoggerEvent evt, const char *message) = 0;
};

template <typename os>
class StreamLogger : public Logger {
private:
    os &stream;
    const constexpr static std::size_t buffer_len = 80U;
    static char buffer[buffer_len];

    void get_time(char *buffer, size_t length)
    {
        std::time_t rawtime;
        std::time(&rawtime);
        auto timeinfo = std::localtime(&rawtime);
        std::strftime(buffer, length, "%Y-%m-%d-%H-%M-%S", timeinfo);
    }

public:
    StreamLogger(os &stream) : stream(stream) { }
    void operator()(LoggerEvent evt, const char *message) override
    {
        get_time(StreamLogger::buffer, StreamLogger::buffer_len);
        stream << "[" << StreamLogger::buffer << " - " << evt << "] "
               << message
               << std::endl;
    }
};

class EmptyLogger : public Logger {
public:
    void operator()(LoggerEvent evt, const char *message) override { }
};

class EventLogger {
private:
    Logger *log;
    LoggerEvent evt;
public:
    EventLogger(Logger *log, LoggerEvent evt)
            : log(log), evt(evt)
    { }

    void operator()(const char *message)
    {
        (*log)(evt, message);
    }
};

template <typename T>
EventLogger &operator<<(EventLogger &logger, const T &t)
{
    std::stringstream s;
    s << t;
    logger(s.str().c_str());
    return logger;
}

#endif //BROTLI_LOGGERS_H
