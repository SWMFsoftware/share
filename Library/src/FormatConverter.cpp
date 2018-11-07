#include "PlotFileIO.h"

int main(int argc, char *argv[]) {

  std::vector<std::string> commandLine;
  for (int i = 0; i < argc; ++i) {
    commandLine.push_back((std::string)(argv[i]));
  }

  FileAgent agent;

  int nFile = commandLine.size() - 1;
  for (int i = 1; i <= nFile; ++i) {
    std::cout << "Processing " << commandLine[i] << "\t\t\t " << i << " / "
              << nFile << std::endl;
    ;
    agent.read(commandLine[i]);
    agent.write();
    // agent.print();
  }
}
