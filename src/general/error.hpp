#ifndef QUEST_ERROR_HPP
#define QUEST_ERROR_HPP

#include <iomanip>
#include <iostream>
#include <sstream>

namespace QUEST {
  
template <typename _T> void Warning(_T msg) {
  if (msg) {
    std::cout << "\n\n" << msg << std::endl;
  }
}

template <typename _T> void Error(_T msg) {
  if (msg) {
    std::cerr << "\n\n" << msg << "\n";
  }
  abort();
}

#define _QUEST_FUNC_NAME __PRETTY_FUNCTION__

#define QUEST_LOCATION \
   "\n ... in function: " << _QUEST_FUNC_NAME << \
   "\n ... in file: " << __FILE__ << ':' << __LINE__ << '\n'

#define QUEST_MESSAGE(MSG, FN)                                   \
  {                                                             \
    std::ostringstream msg_stream;                              \
    msg_stream << std::setprecision(16);                        \
    msg_stream << std::setiosflags(std::ios_base::scientific);  \
    msg_stream << MSG << QUEST_LOCATION;                         \
    QUEST::FN(msg_stream.str().c_str());                         \
  }

#define QUEST_ABORT(MSG) QUEST_MESSAGE("QUEST abort: " << MSG, Error)

#define QUEST_ERROR(MSG) QUEST_MESSAGE("Error: " << MSG, Error)

#define QUEST_WARNING(MSG) QUEST_MESSAGE("Warning: " << MSG, Warning)

#define QUEST_THROW(X, EXCEPTION) \
  if ((X)) {                     \
    throw EXCEPTION;             \
  }

#define QUEST_VERIFY(X, MSG)                                                     \
  if (!(X)) {                                                                   \
    QUEST_MESSAGE("Verification failed: (" << (X) << ") is false\n --> " << MSG, \
                Error);                                                         \
  }

#define QUEST_ASSERT(x, msg)

#define QUEST_SHOW_VAR(var) #var

} // namespace QUEST

#endif  // QUEST_ERROR_HPP
