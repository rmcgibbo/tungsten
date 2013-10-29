#include <string>
#include <stdexcept>
#include <sstream>

namespace Tungsten {

void throwException(const char* file, int line, const std::string& details) {
    std::string fn(file);
    std::string::size_type pos = fn.find_last_of("/\\");
    if (pos+1>=fn.size())
        pos=0;
    std::string filename(fn,(int)(pos+1),(int)(fn.size()-(pos+1)));
    std::stringstream message;
    message << "Assertion failure at "<<filename<<":"<<line << " ";
    if (details.size() > 0)
        message << ".  "<<details;
    throw(std::runtime_error(message.str()));
}

#define ASSERT(cond) {if (!(cond)) throwException(__FILE__, __LINE__, "");};

#define ASSERT(cond) {if (!(cond)) throwException(__FILE__, __LINE__, "");};

#define ASSERT_EQUAL(expected, found) {if (!((expected) == (found))) {std::stringstream details; details << "Expected "<<(expected)<<", found "<<(found); throwException(__FILE__, __LINE__, details.str());}};
 
#define ASSERT_EQUAL_TOL(expected, found, tol) {double _scale_ = std::abs(expected) > 1.0 ? std::abs(expected) : 1.0; if (!(std::abs((expected)-(found))/_scale_ <= (tol))) {std::stringstream details; details << "Expected "<<(expected)<<", found "<<(found); throwException(__FILE__, __LINE__, details.str());}};

#define ASSERT_EQUAL_VEC(expected, found, tol) {double _norm_ = std::sqrt((expected).dot(expected)); double _scale_ = _norm_ > 1.0 ? _norm_ : 1.0; if ((std::abs(((expected)[0])-((found)[0]))/_scale_ > (tol)) || (std::abs(((expected)[1])-((found)[1]))/_scale_ > (tol)) || (std::abs(((expected)[2])-((found)[2]))/_scale_ > (tol))) {std::stringstream details; details << " Expected "<<(expected)<<", found "<<(found); throwException(__FILE__, __LINE__, details.str());}};

#define ASSERT_USUALLY_TRUE(cond) {if (!(cond)) throwException(__FILE__, __LINE__, "(This test is stochastic and may occasionally fail)");};

#define ASSERT_USUALLY_EQUAL_TOL(expected, found, tol) {double _scale_ = std::abs(expected) > 1.0 ? std::abs(expected) : 1.0; if (!(std::abs((expected)-(found))/_scale_ <= (tol))) {std::stringstream details; details << "Expected "<<(expected)<<", found "<<(found)<<" (This test is stochastic and may occasionally fail)"; throwException(__FILE__, __LINE__, details.str());}};

#define ASSERT_VALID_INDEX(index, vector) {if (index < 0 || index >= (int) vector.size()) throwException(__FILE__, __LINE__, "Index out of range");};

}  // namespace Tungsten
