#ifndef CODE_EXPRESSIONML_H
#define CODE_EXPRESSIONML_H

#include <sstream>
#include "geom.h"

namespace ExpressionML {
    class Writer {
        std::ostream& os;

    public:
        explicit Writer(std::ostream& o) : os(o) {}

        void open() {
            os << "<?xml version='1.0'?>"
                  "<!DOCTYPE Expression SYSTEM 'http://www.wolfram.com/XML/notebookml1.dtd'>"
                  "<Expression xmlns:mathematica='http://www.wolfram.com/XML/' xmlns='http://www.wolfram.com/XML/'>";
        }

        void close() {
            os << "</Expression>";
        }

        void openFunction(const std::string& symbol) {
            os << "<Function>";
            writeSymbol(symbol);
        }

        void closeFunction() {
            os << "</Function>";
        }

        void openRule(const std::string& key) {
            openFunction("Rule");
            writeString(key);
        }

        void closeRule() {
            closeFunction();
        }

        void writeNumber(double n) {
            os << "<Number>" << n << "</Number>";
        }

        void writeString(const std::string& string) {
            os << "<String>" << string << "</String>";
        }

        void writeSymbol(const std::string& symbol) {
            os << "<Symbol>" << symbol << "</Symbol>";
        }

        void writePoint(const Point& point) {
            openFunction("List");
            writeNumber(point.x);
            writeNumber(point.y);
            closeFunction();
        }

        void writePoints(const Points& points) {
            openFunction("List");
            for (const auto& point : points)
                writePoint(point);
            closeFunction();
        }

        void writeLine(const Line& line) {
            openFunction("List");
            writePoint(line.origin);
            writePoint(line.direction);
            closeFunction();
        }
    };
}
#endif //CODE_EXPRESSIONML_H
