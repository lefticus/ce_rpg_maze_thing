#include <cstdint>
#include <cstddef>
#include <limits>
#include <utility>
#include <iostream>
#include <string>
#include <array>

/// constexpr random number generator
struct PCG {
  struct pcg32_random_t { std::uint64_t state=0;  std::uint64_t inc=0; };
  pcg32_random_t rng;
  typedef std::uint32_t result_type;

  constexpr result_type operator()() noexcept
  {
    return pcg32_random_r();
  }

  static result_type constexpr min() noexcept
  {
    return std::numeric_limits<result_type>::min();
  }

  static result_type constexpr max() noexcept
  {
    return std::numeric_limits<result_type>::min();
  }

  private:
  constexpr std::uint32_t pcg32_random_r() noexcept
  {
    std::uint64_t oldstate = rng.state;
    // Advance internal state
    rng.state = oldstate * 6364136223846793005ULL + (rng.inc|1);
    // Calculate output function (XSH RR), uses old state for max ILP
    std::uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    std::uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
  }

};

template<typename T, typename Gen>
constexpr auto distribution(Gen &g, T min, T max) noexcept
{
  const auto range = max - min + 1;
  const auto bias_remainder = std::numeric_limits<T>::max() % range;
  const auto unbiased_max = std::numeric_limits<T>::max() - bias_remainder - 1;

  auto r = g();
  for (; r > unbiased_max; r = g());

  return (r % range) + min;
}




//// constexpr map / maze
struct Cell
{
  bool left_open = false;
  bool right_open = false;
  bool up_open = false;
  bool down_open = false;
  bool visited = false;
};

template<typename Data, std::size_t Cols, std::size_t Rows>
struct Array2D
{
  constexpr static auto rows() noexcept { return Rows; }
  constexpr static auto cols() noexcept { return Cols; }

  constexpr Data &operator()(const std::size_t col, const std::size_t row) noexcept
  {
    return m_data[col + (row * Cols)];
  }

  constexpr const Data &operator()(const std::size_t col, const std::size_t row) const noexcept
  {
    return m_data[col + (row * Cols)];
  }

  Data m_data[Cols * Rows];
};

enum class WallType
{
  Empty,
  UpperLeft,
  Vertical,
  Horizontal,
  UpperRight,
  LowerLeft,
  LowerRight,
  RightTee,
  LeftTee,
  UpTee,
  DownTee,
  FourWay,
  Up,
  Down,
  Left,
  Right,
  Visited,
  Dead,
  End,
  Used,
  Monster,
  HiddenMonster
};

template<typename T, std::size_t Size>
struct Stack
{
  T m_data[Size]{};
  std::size_t pos = 0;

  template<typename ... Arg>
    void constexpr emplace_back(Arg &&... arg)
    {
      m_data[pos] = T{std::forward<Arg>(arg)...};
      ++pos;
    }
  
  constexpr T pop_back()
  {
    --pos;
    return m_data[pos];
  }

  constexpr bool empty() const
  {
    return pos == 0;
  }

  constexpr std::size_t size() const
  {
    return pos;
  }
};

struct Loc
{
  std::size_t col=0;
  std::size_t row=0;

  constexpr Loc(const std::size_t t_col, const std::size_t t_row)
    : col(t_col), row(t_row)
  {
  }

  constexpr Loc() = default;
};

template<std::size_t num_cols, std::size_t num_rows>
constexpr Array2D<Cell, num_cols, num_rows> make_maze(const std::uint64_t seed)
{

  PCG pcg{0, seed};

  Array2D<Cell, num_cols, num_rows> M;



  // Set starting row and column
  std::size_t c = 0;
  std::size_t r = 0;
  Stack<Loc, num_cols * num_rows> history;
  history.emplace_back(c,r); // The history is the stack of visited locations

  // Trace a path though the cells of the maze and open walls along the path.
  // We do this with a while loop, repeating the loop until there is no history, 
  // which would mean we backtracked to the initial start.
  while (!history.empty()) 
  {
    M(c, r).visited = true;

    Stack<char, 4> check{};

    // check if the adjacent cells are valid for moving to
    if (c > 0 && M(c-1, r).visited == false) {
      check.emplace_back('L');
    }
    if (r > 0 && M(c, r-1).visited == false) {
      check.emplace_back('U');
    }
    if (c < num_cols-1 && M(c+1, r).visited == false) {
      check.emplace_back('R');
    }
    if (r < num_rows-1 && M(c, r+1).visited == false) {
      check.emplace_back('D');
    }

    if (!check.empty()) { // If there is a valid cell to move to.
      // Mark the walls between cells as open if we move

      history.emplace_back(c,r);

      for (auto num_pops = distribution(pcg, std::size_t(0), check.size() - 1); num_pops > 0; --num_pops)
      {
        check.pop_back();
      }

      switch (check.pop_back()) {
        case 'L':
          M(c, r).left_open = true;
          --c;
          M(c, r).right_open = true;
          break;

        case 'U':
          M(c, r).up_open = true;
          --r;
          M(c, r).down_open = true;
          break;

        case 'R':
          M(c, r).right_open = true;
          ++c;
          M(c, r).left_open = true;
          break;

        case 'D':
          M(c, r).down_open = true;
          ++r;
          M(c, r).up_open = true;
          break;
      }
    } else {
      // If there are no valid cells to move to.
      // retrace one step back in history if no move is possible
      const auto last = history.pop_back();
      c = last.col;
      r = last.row;
    }
  }

  // Open the walls at the start and finish
  M(0,0).left_open = true;
  M(num_cols-1, num_rows-1).right_open = true;

  return M;
}



template<std::size_t num_cols, std::size_t num_rows>
constexpr Array2D<WallType, num_cols*2+1, num_rows*2+1> render_maze(const Array2D<Cell, num_cols, num_rows> &maze_data)
{
  Array2D<WallType, num_cols*2+1, num_rows*2+1> result{};


  for (std::size_t col = 0; col < num_cols; ++col)
  {
    for (std::size_t row = 0; row < num_rows; ++row)
    {
      const auto render_col = col * 2 + 1;
      const auto render_row = row * 2 + 1;

      const auto &cell = maze_data(col, row);

      // upper
      if (!cell.up_open) { result(render_col,render_row-1) = WallType::Horizontal; }

      // left
      if (!cell.left_open) { result(render_col-1,render_row) = WallType::Vertical; }

      // right
      if (!cell.right_open) { result(render_col+1,render_row) = WallType::Vertical; }

      // lower
      if (!cell.down_open) { result(render_col,render_row+1) = WallType::Horizontal; }
    }
  }

  for (std::size_t col = 0; col < result.cols(); col += 2)
  {
    for (std::size_t row = 0; row < result.rows(); row += 2)
    {
      const auto &cell = result(col, row);

      const auto up     = (row == 0)?false:result(col, row-1) != WallType::Empty;
      const auto left   = (col == 0)?false:result(col-1, row) != WallType::Empty;
      const auto right  = (col == num_cols * 2)?false:result(col+1, row) != WallType::Empty;
      const auto down   = (row == num_rows * 2)?false:result(col, row+1) != WallType::Empty;

      if (up && right && down && left) {
        result(col, row) = WallType::FourWay;
      }
      if (up && right && down && !left) {
        result(col, row) = WallType::RightTee;
      }
      if (up && right && !down && left) {
        result(col, row) = WallType::UpTee;
      }
      if (up && !right && down && left) {
        result(col, row) = WallType::LeftTee;
      }
      if (!up && right && down && left) {
        result(col, row) = WallType::DownTee;
      }

      if (up && right && !down && !left) {
        result(col, row) = WallType::LowerLeft;
      }
      if (up && !right && !down && left) {
        result(col, row) = WallType::LowerRight;
      }
      if (!up && !right && down && left) {
        result(col, row) = WallType::UpperRight;
      }
      if (!up && right && down && !left) {
        result(col, row) = WallType::UpperLeft;
      }
      if (!up && right && !down && left) {
        result(col, row) = WallType::Horizontal;
      }
      if (up && !right && down && !left) {
        result(col, row) = WallType::Vertical;
      }


      if (up && !right && !down && !left) {
        result(col, row) = WallType::Up;
      }
      if (!up && right && !down && !left) {
        result(col, row) = WallType::Right;
      }
      if (!up && !right && down && !left) {
        result(col, row) = WallType::Down;
      }
      if (!up && !right && !down && left) {
        result(col, row) = WallType::Left;
      }
    }
  }


  return result;
}

enum struct Command {
  Left,
  Right,
  Up,
  Down,
  Stay
};

template<typename T, std::size_t CommandSize>
constexpr auto execute_commands(T maze, std::array<Command, CommandSize> commands)
{
    std::size_t col = 0;
    std::size_t row = 1;
    for (const auto &command : commands) {
      bool fail = false;
      switch (command) {
          case Command::Left:
            if (col == 0) {
                fail = true;
            } else {
                --col;
            }
            break;
          case Command::Right:
            ++col;
            break;
          case Command::Up:
            --row;
            break;
          case Command::Down:
            ++row;
            break;
          case Command::Stay:
            break;
      }
      if (maze(col, row) == WallType::HiddenMonster) {
        maze(col, row) = WallType::Monster;
        break;
      } else if (maze(col, row) != WallType::Empty || fail) {
        maze(col, row) = WallType::Dead;
        break;
      } else if (col == maze.cols() - 1 && row == maze.rows() - 2) {
        maze(col, row) = WallType::End;
      } else {
        maze(col, row) = WallType::Used;
      }

    }
    return maze;
}

template<typename Maze>
constexpr auto add_monsters(Maze maze, const std::uint64_t seed)
{
    PCG pcg{0, seed};

    for (std::size_t row = 1; row < maze.rows() - 1; ++row) {
        for (std::size_t col = 1; col < maze.cols() - 1; ++col) {
            if (const auto type = maze(col, row); 
                type == WallType::Vertical || type == WallType::Horizontal ) 
            {
                if (distribution(pcg, 1, 20) == 2) {
                    maze(col, row) = WallType::HiddenMonster;
                }
            }
        }
    }

    return maze;
}

template<typename T>
auto print_maze(const T &maze)
{

  for (std::size_t row = 0; row < maze.rows(); ++row)
  {
    for (std::size_t col = 0; col < maze.cols(); ++col)
    {
      const auto square = [&](){
        switch (maze(col, row)) {
          case WallType::Empty:         return " ";
          case WallType::UpperLeft:     return "┌";
          case WallType::Vertical:      return "│";
          case WallType::Horizontal:    return "─";
          case WallType::UpperRight:    return "┐";
          case WallType::LowerLeft:     return "└";
          case WallType::LowerRight:    return "┘";
          case WallType::RightTee:      return "├";
          case WallType::LeftTee:       return "┤";
          case WallType::UpTee:         return "┴";
          case WallType::DownTee:       return "┬";
          case WallType::FourWay:       return "┼";
          case WallType::Up:            return "│";
          case WallType::Down:          return "│";
          case WallType::Left:          return "─";
          case WallType::Right:         return "─";
          case WallType::Visited:       return "·";
          case WallType::Used:          return "*";
          case WallType::Dead:          return "🕱";
          case WallType::End:           return "😄";
          case WallType::Monster:       return "👾";
          case WallType::HiddenMonster: return " ";
        }
        throw "Unhandled";
      }();

      std::cout << square;
    }
    std::cout << '\n';
  }

}
