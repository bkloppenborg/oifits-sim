/*
 * Read.cpp
 *
 *  Created on: Feb 18, 2011
 *      Author: bkloppenborg
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <algorithm>

#include "ReadTextFile.h"

using namespace std;

// A simple string splitting function for use with delimiters
void StringSplit(string str, string delim, vector<string> & results)
{
    int cutAt;
    while( (cutAt = str.find_first_of(delim)) != str.npos )
    {
        if(cutAt > 0)
            results.push_back(str.substr(0,cutAt));

        str = str.substr(cutAt+1);
    }

    // Push any remaining characters on to the back of the vector.
    if(str.length() > 0)
        results.push_back(str);
}

string StripWhitespace(string str)
{
    str.erase(remove_if(str.begin(), str.end(), ::isspace), str.end());
    return str;
}


void StripWhitespace(vector<string> & strings)
{
    // Remove any remaining white space.
    for(unsigned int i = 0; i < strings.size(); i++)
    {
        strings[i] = StripWhitespace(strings[i]);
    }
}

/// Reads in a file, returns non-comment lines as a vector of strings
/// If the file cannot be opened, an exception is thrown with the message contained in error_message
vector<string> ReadFile(string filename, string comment_chars, string error_message)
{
    ifstream infile;
    infile.open(filename.c_str());
    vector < string > lines;

    if (infile.is_open())
    {
        string line;

        while (!infile.eof())
        {
            getline(infile, line);
            while ((line.size() == 0 || comment_chars.find(line[0]) != string::npos)
                   && !infile.eof())
            {
                getline(infile, line);
            }
            if (!infile.eof())
                lines.push_back(line);
        }
        infile.close();
    }
    else
    {
        /// \exception runtime_error Error opening array file
        /// The array file could not be located.  It is likely that the user just specified an invalid path.
        throw std::runtime_error(error_message);
    }

    return lines;
}

// Tokenizes a line.
vector<string> Tokenize(string line)
{
    istringstream lineStream(line);

    vector < string > tokens;
    while (lineStream)
    {
        string item;

        lineStream >> item;
        tokens.push_back(item);
    }

    return tokens;
}

// See if a file exists.  Returns true if it does.
bool FileExists(string filename)
{
  ifstream ifile(filename.c_str());
  return bool(ifile);
}
