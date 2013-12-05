#ifndef _ASTRALCAT_LOGGER_
#define _ASTRALCAT_LOGGER_

#include <iostream>
#include <memory>
#include <boost/format.hpp>
#include <boost/utility.hpp>
#include <boost/lexical_cast.hpp>
#include <stdexcept>


namespace astralcat {

    class Logger {
    public:
        enum level_t {DEBUG, INFO, WARN, ERROR, FATAL};
        Logger(std::ostream &os = std::cerr, level_t level = INFO);
        template <typename... Args> Logger &debug (const char *format, const Args &... args) { return log(DEBUG, format, args...); }
        template <typename... Args> Logger &info  (const char *format, const Args &... args) { return log(INFO,  format, args...); }
        template <typename... Args> Logger &warn  (const char *format, const Args &... args) { return log(WARN,  format, args...); }
        template <typename... Args> Logger &error (const char *format, const Args &... args) { return log(ERROR, format, args...); }
        template <typename... Args> Logger &fatal (const char *format, const Args &... args) { return log(FATAL, format, args...); }

        friend class Indent;

        class Indent : boost::noncopyable {
            Logger *logger;
        public:
            Indent(Logger *logger) : logger(logger) {
                logger->retain();
            }
            Indent(const Indent &indent) {
                logger->retain();
            }
            ~Indent() { logger->release(); }
        };

        Indent indent() {
            return Indent(this);
        }

        operator Indent() { return indent(); }

    private:
        struct Impl;
        std::shared_ptr<Impl> pimpl;

        boost::format &apply_args(boost::format &fmt) {
            return fmt;
        }

        template <typename T, typename... Args>
        boost::format &apply_args(boost::format &fmt, const T &a0, const Args &... args) {
            return apply_args(fmt % a0, args...);
        }

        template <typename... Args>
        Logger &log(level_t level, const char *format, const Args &... args) try {
            boost::format fmt(format);
            apply_args(fmt, args...);
            out(level, fmt);
            return *this;
        }
        catch (const std::exception &e) {
            std::cerr << "logger error: " << e.what() << std::endl;
            return *this;
        }

        void out(level_t level, boost::format &fmt);
        void retain();
        void release();

    };

    extern Logger logger;

}

#endif
