#include "Simulator.h"

double sinc( double number )
{
  if (abs(number) < 1.0e-10)
    return 1.0;
  return sin(number) / number;
}

// A simple string splitting function for use with delimiters
void StringSplit(string str, string delim, vector<string> results)
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

void StripWhitespace(string str)
{
    str.erase(remove_if(str.begin(), str.end(), isspace), str.end());
    return str;
}

void StripWhitespace(vector<string> strings)
{
    // Remove any remaining white space.
    for(int i = 0; i < strings.size(); i++)
    {
        strings[i] = StripWhitespace(strings[i]);
    }
}
