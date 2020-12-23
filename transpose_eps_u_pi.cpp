#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <string>

using namespace std;

int main(int argc, char *argv[])
{
	if (argc < 2) exit(8);

	// read path to input file from command line
	string path_to_file = string(argv[1]);

	int xsize = -1, ysize = -1;

	vector<string> all_lines;
	
	// then read in file itself
	ifstream infile(path_to_file.c_str());
	if (infile.is_open())
	{
		string line;
		int lineCount = 0;
		while ( getline (infile, line) )
		{
			// skip header line; just process and print it back out
			if ( lineCount++ < 1 )
			{
				string dummy;
				istringstream issHeader(line);
				issHeader >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy
                          >> xsize >> dummy >> ysize >> dummy >> dummy >> dummy
                          >> dummy >> dummy >> dummy;
				cout << line << endl;
				continue;
			}
			// skip empty lines; don't print them back out since they're not
			// included in IP-Glasma default output
			if ( line.empty() ) continue;

			// otherwise, if the line contains stuff, store it
			all_lines.push_back( line );
		}
	}

	// check for problems
	if ( xsize < 1 or ysize < 1 or all_lines.size() != xsize*ysize )
	{
		cerr << "Problem with reading in file!  Aborting..." << endl;
		exit(8);
	}
	else
	{
		// print 
		for (int ix = 0; ix < xsize; ix++)
		for (int iy = 0; iy < ysize; iy++)
			cout << all_lines[iy*xsize + ix] << endl;
		cout << endl;
	}
	
	infile.close();

	return (0);
}


