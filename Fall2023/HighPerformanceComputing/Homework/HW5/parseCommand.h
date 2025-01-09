#ifndef PARSE_COMMAND_H
#define PARSE_COMMAND_H

#include <string>
using std::string;

// ------------------------------------------------------------------------------------
// Parse a command from a string "line" and return a integer "value"
// Example:
// line= "-Nx=10" (input)
// command= "-Nx=" (input)
// value= 10 (output)
// Return value : 1 if command was found in line
// 0 if not found
// ------------------------------------------------------------------------------------
int parseCommand( const string & line, const string & command, int & value) {
    int len =command.length();
    if( line.substr(0,len)==command ) {
        sscanf(line.substr(len).c_str(),"%d",&value);
        printf("parseCommand: SETTING %s%d\n",command.c_str(),value);
        return 1;
    }
    return 0;
}

// Parse a double
int parseCommand( const string & line, const string & command, double & value) {
    int len =command.length();
    if( line.substr(0,len)==command ) {
        sscanf(line.substr(len).c_str(),"%le",&value);
        printf("parseCommand: SETTING %s%e\n",command.c_str(),value);
        return 1;
    }
    return 0;
}

// Parse a string
int parseCommand( const string & line, const string & command, string & value) {
    int len =command.length();
    if( line.substr(0,len)==command ) {
        value = line.substr(len);
        printf("parseCommand: SETTING %s%s\n",command.c_str(),value.c_str());
        return 1;
    }
    return 0;
}


#endif
