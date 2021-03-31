#ifndef READPARAM_H
#define READPARAM_H

#include <climits>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

class ReadParam {
public:
  std::stringstream ss;
  bool isVerbose;

  std::string commandSuffix = std::string();

public:
  ReadParam() { isVerbose = true; }
  ~ReadParam() {}
  ReadParam &operator=(const std::string &stringIn) {
    // ss << stringIn;
    ss.str(stringIn);
    return (*this);
  }
  ReadParam &operator=(const ReadParam &other) { return (*this); }

  void set_verbose(bool in) { isVerbose = in; }
  //==========================================

  void set_command_suffix(std::string &in) { commandSuffix = in; }

  inline bool get_next_command(std::string &id) {
    id.clear();
    ss.ignore(INT_MAX, '#');
    while (ss.good()) {
      ss.unget();

      // Example of commandline:
      // 1) #COMMAND
      // 2) #COMMAND_GRID0_GRID1

      std::string commandline;
      ss >> commandline;
      std::string::size_type pos1;
      // Check if this command has suffixes
      pos1 = commandline.find_first_of("_");
      id = commandline.substr(0, pos1);

      if (pos1 != std::string::npos) {
        std::string tmp = "_" + commandSuffix;
        if (commandline.find(tmp) == std::string::npos) {
          // Go to next command.
          ss.ignore(INT_MAX, '#');
          continue;
        }
      }
      if (isVerbose)
        std::cout << "\n"
                  << "PC: " << id << " " << commandSuffix << std::endl;
      ss.ignore(INT_MAX, '\n');
      return true;
    }

    return false;
  }
  //==============================

  template <class T> void read_var(std::string description, T &var) {

    if (ss.good())
      ss >> var; // If ss is empty, ss.good() is still true.

    if (ss.good()) {
      // If ss is empty, ss.good() is false now.
      // If ss is the last line, ss.good() is still true is there is a
      // '\n' at the end of each line.
      ss.ignore(INT_MAX, '\n');
      if (isVerbose) {
        std::cout << "PC: " << std::left << std::setw(40) << var << description
                  << std::endl;
        return;
      }
    } else {
      if (isVerbose)
        std::cout << "Type T: can not find input parameter. Abort!"
                  << std::endl;
      std::abort();
    }
  }
  //=======================================================

  inline void read_var(std::string description, bool &var) {
    std::string text;

    if (ss.good())
      ss >> text;

    if (ss.good()) {
      ss.ignore(INT_MAX, '\n');
      var = false;
      if (text == "T")
        var = true;
      if (isVerbose) {
        std::cout << "PC: " << std::left << std::setw(40) << text << description
                  << std::endl;
        return;
      }
    } else {
      if (isVerbose)
        std::cout << "Bool: can not find input parameter. Abort!" << std::endl;
      std::abort();
    }
  }
  //==========================================================

  inline void read_var(std::string description, std::string &var) {
    // Read input until \t or next line or three successive white spaces.
    std::string text;

    std::string temp;
    bool isEmpty = false;
    if (ss.good()) {
      std::getline(ss, temp);
      std::string delimiter = "\t\n ";
      // Make sure temp is not empty or does not contain white characters only.
      isEmpty = temp.empty() ||
                temp.find_first_not_of(delimiter) == std::string::npos;
    }

    if (!isEmpty) {
      std::string::size_type pos;

      // Replace the first three successive white spaces with tab.
      pos = temp.find("   "); // Three white spaces.
      if (pos != std::string::npos) {
        temp.replace(pos, pos + 3, "\t");
      }

      std::string delimiter = "\t\n"; // Not include space ' '.
      pos = temp.find_first_of(delimiter);
      if (pos != std::string::npos)
        temp.erase(pos);
      // Do not have '\'. Exclude unvisible junk characters.
      std::string char0 =
          "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz01"
          "23456789_=+-.~&*[]|{}@!#$%^()/?<>,;:";
      var = temp.substr(0, temp.find_last_of(char0) + 1);
      if (isVerbose) {
        std::cout << "PC: " << std::left << std::setw(40) << var << description
                  << std::endl;
        return;
      }
    } else {
      if (isVerbose)
        std::cout << "String: can not find input parameter. Abort!"
                  << std::endl;
      std::abort();
    }
  }
  //=========================================

  inline void skip_lines(int nlines) {
    for (int i = 0; i < nlines; i++) {
      ss.ignore(INT_MAX, '\n');
    }
  }
};

extern bool isProc0;

inline void char_to_string(std::string &ss, char *chararray, int length,
                           int linelength) {
  for (int i = linelength; i < length; i += linelength)
    chararray[i - 1] = '\n';
  ss.assign(chararray, length);
}

inline std::string char_to_string(char *chars, int length, int linelength) {
  for (int i = linelength; i < length; i += linelength)
    chars[i - 1] = '\n';
  return std::string(chars, length);
}

inline std::stringstream *char_to_stringstream(char *chararray, int length,
                                               int linelength, int iProc) {

  isProc0 = iProc == 0;

  std::stringstream *ss;
  for (int i = linelength; i < length; i += linelength)
    chararray[i - 1] = '\n';
  ss = new std::stringstream;
  ss->write(chararray, length);
  return (ss);
}

template <class T>
void read_var(std::stringstream *ss, std::string description, T *var) {
  if (*ss) {
    *ss >> *var;
    ss->ignore(INT_MAX, '\n');
    if (isProc0) {
      std::cout << "PC: " << std::left << std::setw(40) << *var << description
                << std::endl;
      return;
    }
  } else {
    std::cout << " Can not find input parameter. Abort!" << std::endl;
    std::abort();
  }
}

inline void read_var(std::stringstream *ss, std::string description,
                     bool *var) {
  std::string text;
  if (*ss) {
    *ss >> text;
    ss->ignore(INT_MAX, '\n');
    *var = false;
    if (text == "T")
      *var = true;
    if (isProc0) {
      std::cout << "PC: " << std::left << std::setw(40) << text << description
                << std::endl;
      return;
    }
  } else {
    std::cout << " Can not find input parameter. Abort!" << std::endl;
    std::abort();
  }
}

inline void read_var(std::stringstream *ss, std::string description,
                     std::string *var) {
  // Read input until \t or next line or three successive white spaces.
  std::string text;
  if (*ss) {
    std::string temp;
    std::getline(*ss, temp);
    std::string::size_type pos;

    // Replace the first three successive white spaces with tab.
    pos = temp.find("   "); // Three white spaces.
    if (pos != std::string::npos) {
      temp.replace(pos, pos + 3, "\t");
    }

    std::string delimiter = "\t\n"; // Not include space ' '.
    pos = temp.find_first_of(delimiter);
    if (pos != std::string::npos)
      temp.erase(pos);
    // Do not have '\'. Exclude unvisible junk characters.
    std::string char0 = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz01"
                        "23456789_=+-.~&*[]|{}@!#$%^()/?<>,;:";
    (*var) = temp.substr(0, temp.find_last_of(char0) + 1);
    if (isProc0) {
      std::cout << "PC: " << std::left << std::setw(40) << (*var) << description
                << std::endl;
      return;
    }
  } else {
    std::cout << " Can not find input parameter. Abort!" << std::endl;
    std::abort();
  }
}

inline void get_next_command(std::stringstream *ss, std::string *id) {
  *id = "";
  ss->ignore(INT_MAX, '#');
  ss->unget();
  *ss >> *id;
  if (isProc0 && id->find('#') != std::string::npos)
    std::cout << "\n"
              << "PC: " << *id << std::endl;
  ss->ignore(INT_MAX, '\n');
}

inline void skip_lines(std::stringstream *ss, int nlines) {
  for (int i = 0; i < nlines; i++) {
    ss->ignore(INT_MAX, '\n');
  }
}

#endif
