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

        void openFunction(std::string const& symbol) {
            os << "<Function>";
            writeSymbol(symbol);
        }

        void closeFunction() {
            os << "</Function>";
        }

        void openRule(std::string const& key) {
            openFunction("Rule");
            writeString(key);
        }

        void closeRule() {
            closeFunction();
        }

        void writeNumber(double n) {
            os << "<Number>" << n << "</Number>";
        }

        void writeString(std::string const& string) {
            os << "<String>" << string << "</String>";
        }

        void writeSymbol(std::string const& symbol) {
            os << "<Symbol>" << symbol << "</Symbol>";
        }

        void writePoint(Point const& point) {
            openFunction("List");
            writeNumber(point.x);
            writeNumber(point.y);
            closeFunction();
        }

        void writePoints(Points const& points) {
            openFunction("List");
            for (auto const& point : points)
                writePoint(point);
            closeFunction();
        }

        void writeLine(Line const& line) {
            openFunction("List");
            writePoint(line.origin);
            writePoint(line.direction);
            closeFunction();
        }
    };
}
#endif //CODE_EXPRESSIONML_H
