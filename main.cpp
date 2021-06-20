/*
    Christian Valencia
    2275944
    valen193@mail.chapman.edu
    CPSC 350 - Data Structures
    #Assignment 1
*/

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <random>
#include <time.h>

using namespace std;

///Normalizes a string to be ALL CAPS
string UpperString(string str)
{
    string up = str;
    int i = 0;
    char c;
    while (up[i])
    {
      up[i] = toupper(up[i]);
      i++;
    }
    return up;
}

//returns true if the string consists only of A,T,C,G
bool validString(string str)
{
  char c;
  for (int i = 0; i < (str.size()-1); i++)
  {
    c = toupper(str[i]);
    if (c != 'A')
    {
      if (c != 'C')
      {
        if (c != 'G')
        {
          if (c != 'T')
          {
            return 0;
          }
        }
      }
    }
  }
  return 1;
}

//A simple count of a file's valid lines
int countLines(ifstream& file)
{
    int num = 0;
    string line;
    if(file.is_open())
    {
      //These two files reset the getline function
      file.clear();
      file.seekg(0);

      while (getline (file, line))
      {
        if(validString(line))
        {
          num ++;
        }
      }
    }
    return num;
}

//Finds the sum of the lenghts of a file's lines
int findSum(ifstream& file)
{
    int s = 0;
    string line;
    if(file.is_open())
    {
      file.clear();
      file.seekg(0);

      while (getline (file, line))
      {
        if(validString(line))
        {
          s += (line.length() - 1);
        }
      }
    }
    return s;
}

//Quick maths
float findMean(int sum, int num)
{
   return sum/num;
}

//Taking the difference of the mean squared of each value
//and divides by the number of lines.
double findVariance(ifstream& file, int num, float mean)
{
  double var = 0;
  string line;
  if(file.is_open())
  {
    file.clear();
    file.seekg(0);

    while (getline (file, line))
    {
      if(validString(line))
      {
        var += pow(((line.size() - 1) - mean), 2.0);
      }
    }
  }
  return var/num;
}

//Standard deviation
double findSD(double var)
{
    return sqrt(var);
}

//Enter a captial nucleotide to find the probability
float singleProb(ifstream& file, char nucleotide, ofstream& out)
{
  float total = 0;
  float matches = 0;
  string line;
  char c;

  if(file.is_open())
  {
    file.clear();
    file.seekg(0);

    while (getline (file, line))
    {
      if(validString(line))
      {
        for (int i = 0; i < (line.size()-1); i++)
        {
          //Here we keep track of the total and the number of matches
          total++;
          c = toupper(line[i]);
          if (c == nucleotide) matches++;
        }
      }
    }
  }
  //Then calculate the probability
  float prob = matches/total * 100;

//If there is a file open we print directly out
//to it and return the value for later use
  if(out.is_open())
  {
    out << "There is a: " << prob << "% chance of " << nucleotide << endl;
  }

  return prob;
}

//Like the previous function but takes in two nucleotides
void biGramProb(ifstream& file, string bi, ofstream& out)
{
  float total = 0;
  float matches = 0;
  string line;
  string sub;

  if(file.is_open())
  {
    file.clear();
    file.seekg(0);

    while (getline (file, line))
    {
      if(validString(line))
      {
        line = UpperString(line);
        //Divides up string
        for (int i = 0; i < (line.size()-1); i++)
        {
          sub = line.substr(i,2);
          if (sub == bi)
          {
            matches++;
          }
          total ++;
          i++;
        }
      }
    }
  }
  float prob = matches/total * 100;
  if(out.is_open())
  {
    out << "There is a: " << prob << "% chance of bi-gram " << bi << endl;
  }
}

//Itterate through all possible bigrams with some string manipulation.
void findAllBiGrams(ifstream& file, ofstream& out)
{
  int i,j;
  string bi;
  float prob;
  for(i = 0; i < 4; i++){
    for(j = 0; j < 4; j++){
      if (i == 0) bi.replace(0,0,"A");
      if (i == 1) bi.replace(0,0,"T");
      if (i == 2) bi.replace(0,0,"C");
      if (i == 3) bi.replace(0,0,"G");

      if (j == 0) bi.replace(1,1,"A");
      if (j == 1) bi.replace(1,1,"T");
      if (j == 2) bi.replace(1,1,"C");
      if (j == 3) bi.replace(1,1,"G");

      biGramProb(file, bi, out);
      //reset the string
      bi = "";
    }
  }

}

string probability(double num, float aProb, float tProb, float cProb, float gProb)
{
  double r = num;
  if ((r -= aProb) < 0) return "A";
  if ((r -= tProb) < 0) return "T";
  if ((r -= cProb) < 0) return "C";
  if ((r -= gProb) < 0) return "G";
}

//Takes all the previous data and creates its own set of sequences
void generateNewDNA(int amount, float mean, double sd, float aProb, float tProb, float cProb, float gProb, ofstream& out)
{
  if (out.is_open())
  {
    int repeat = amount;
    double a,b,C,D,r;
    double pi = 3.1415926535897;
    string nucleotide;

    //initialize seed for randoms
    srand (time(NULL));
    default_random_engine random(time(0));
    uniform_real_distribution<double> distribution(0.0,1.0);

    while (repeat > 0)
    {
      //create random number (0,1)
      a = distribution(random);
      b = distribution(random);

      //ensure number is [0,1)
      while (a == 0) {
        a = distribution(random);
      }
      while (b == 0) {
        b = distribution(random);
      }


      C = (sqrt(-2 * log(a))* cos(2*pi*b));
      D = sd * C + mean;

      for(int i = 0; i < D; i++)
      {
        r = (rand() % 100 + 1);
        nucleotide = probability(r, aProb, tProb, cProb, gProb);
        out << nucleotide;
      }
      if(D > 0)
      {
        out << endl;
        repeat--;
      }
    }
  }
}

//Takes in a filepath and runs all calculations and outputs to Valencia.out
void analyze(string filepath)
{
  cout << "Reading from " << filepath << endl;
  ofstream fileOut ("results.txt", ios::app);
  if (fileOut.is_open())
  {
    ifstream fileIn (filepath);
    if(fileIn.is_open())
    {

      fileOut << "Christian Valencia" << endl;
      fileOut << "2275944" << endl;
      fileOut << "valen193@mail.chapman.edu" << endl;
      fileOut << "CPSC 350 - Data Structures" << endl;
      fileOut << "#Assignment 1" << endl;
      fileOut << endl;
      fileOut << "Begin analysis of " << filepath << endl;
      fileOut << endl;
      int num = countLines(fileIn);
      fileOut << "The total number of valid lines is: " << num << endl;

      int sum = findSum(fileIn);
      fileOut << "The sum of the lengths is: " << sum << endl;

      float mean = findMean(sum,num);
      fileOut << "The mean of the lengths is: " << mean << endl;

      double variance = findVariance(fileIn,num,mean);
      fileOut << "The variance of the lengths is: " << variance << endl;

      double sd = findSD(variance);
      fileOut << "The standard deviation of the lengths is: " << sd << endl;

      float aProb = singleProb(fileIn, 'A', fileOut);
      float tProb = singleProb(fileIn, 'T', fileOut);
      float cProb = singleProb(fileIn, 'C', fileOut);
      float gProb = singleProb(fileIn, 'G', fileOut);

      findAllBiGrams(fileIn, fileOut);
      fileOut << endl;
      generateNewDNA(1000,mean,sd,aProb,tProb,cProb,gProb,fileOut);

      fileOut << endl;
      fileOut << "End of file" << endl;
      fileOut << endl;

      fileIn.close();
      cout << "Results have be saved to results.txt" << endl;
    }
    else
    {
      cout << "Unrecognized file path. Please try again." << endl;
    }

    fileOut.close();
  }
}

int main(int argc, char** argv)
{
  if((argc > 0) && (argv[1] != NULL))
  {
    bool done = false;
    string io;

    string filepath =  argv[1];

    //Main function
    analyze(filepath);

    //Looping mechanisms
    while (!done)
    {
      cout << "Would you like to analyze another file? [Y/N]" << endl;
      cin >> io;
      io = UpperString(io);
      if (io == "N")
      {
        cout << "Goodbye!" << endl;
        done = true;
        break;
      }
      else if (io == "Y")
      {
        cout << "Please enter the another file path: ";
        cin >> filepath;
        analyze(filepath);
      }
      else
      {
        cout << "Unrecognized input. Please try again." << endl;
      }
    }
  }
  else
  {
      cout << "Please provide the file path for a text document" << endl;
      cout << "USAGE ./main.out <file-path>" << endl;
  }
  cout << "The program exited successfully" << endl;
  return 0;
}
