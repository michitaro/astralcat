#include "Logger.h"
#include <stdlib.h>
#include <string.h>


namespace astralcat {

    const char *const CLEAR  = "\033[0m";
    const char *const BLACK  = "\033[30m";
    const char *const RED    = "\033[31m";
    const char *const GREEN  = "\033[32m";
    const char *const YELLOW = "\033[33m";
    const char *const BLUE   = "\033[34m";
    const char *const PURPLE = "\033[35m";
    const char *const CYAN   = "\033[36m";
    const char *const WHITE  = "\033[37m";


    struct Logger::Impl {
        int               indent;
        Logger::level_t   level;
        std::ostream      &os;

        Impl(std::ostream &os, level_t level) :
            level(level),
            os(os)
        {
        }

        void write_prefix(Logger::level_t level) {
            switch (level) {
                case DEBUG:  os << BLUE   << "[D] ";   break;
                case INFO:   os << GREEN  << "[I] ";   break;
                case WARN:   os << YELLOW << "[W] ";   break;
                case ERROR:  os << RED    << "[E] ";   break;
                case FATAL:  os << RED    << "[F] ";   break;
            }
        }

        void out(Logger::level_t level, boost::format &fmt) {
            if (level < this->level) {
                return;
            }
            write_prefix(level);
            for (int i = 0;  i < indent;  i++) {
                os << "  ";
            }
            os << fmt << CLEAR << std::endl;
        }

    };

    Logger::Logger(std::ostream &os, level_t level) : pimpl(std::make_shared<Logger::Impl>(os, level))
    {
    }

    void Logger::out(Logger::level_t level, boost::format &fmt) {
        pimpl->out(level, fmt);
    }

    void Logger::retain() {
        pimpl->indent++;
    }

    void Logger::release() {
        pimpl->indent--;
    }


    namespace {
        Logger default_logger() {
            Logger::level_t level = Logger::INFO;
            if (const char *env = getenv("ASTRALCAT_LOGLEVEL")) {
                if (strcmp(env, "debug") == 0) level = Logger::DEBUG;
                if (strcmp(env, "info")  == 0) level = Logger::INFO;
                if (strcmp(env, "warn")  == 0) level = Logger::WARN;
                if (strcmp(env, "error") == 0) level = Logger::ERROR;
                if (strcmp(env, "fatal") == 0) level = Logger::FATAL;
            }
            return Logger(std::cerr, level);
        }
    }

    Logger logger = default_logger();;

}
