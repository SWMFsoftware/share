#ifndef READPARAM_H
#define READPARAM_H

#include <iomanip>
#include <iostream>
#include <sstream>
#include <climits>
#include <cstdlib>

extern bool isProc0;

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
                     string *var) {
  // Read input until \t or next line or three successive white spaces.
  std::string text;
  if (*ss) {
    string temp;
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
