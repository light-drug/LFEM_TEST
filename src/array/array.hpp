#ifndef QUEST_ARRAY_HPP
#define QUEST_ARRAY_HPP

#include "config.hpp"
#include <array>

namespace QUEST{

  using Array = std::array<real_t, QUEST_SPACE_DIM>;
  using Coordinate = std::array<real_t, QUEST_SPACE_DIM>;
  using IntArray = std::array<int, QUEST_SPACE_DIM>;

} // namespace QUEST

#endif  // QUEST_ARRAY_HPP 