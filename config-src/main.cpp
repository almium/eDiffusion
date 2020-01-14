#include <iostream>
#include <stdio.h>
using namespace std;

#include "config.h"

int main(int argc, char* argv[], char* envp[])
{
	cout << "config demo program" << endl;
	cout << "(c) 2008 Kai G. Schwebke" << endl;
	cout << "This code is licensed under LGPL and comes with no warranty;" << endl;
	cout << "see COPYING for details." << endl;
	cout << endl;

	// read config file with environment variable expansion support
	Config config("demo.config", envp);
	
	// basic usage of properties
	if (config.pBool("showHello")) {
		int cnt = config.pInt("helloCount");
		string msg = config.pString("helloMessage");
		for (int i = 0; i < cnt; ++i) {
			cout << msg << endl;
		}
		cout << endl;
	}

	
	printf(" pi = %2.6f \n",config.pDouble("pi"));


	// properties with expanded names (no difference in using)
	cout << "tempFolder    = '" << config.pString("tempFolder") << "'" << endl;
	cout << "tempSubFolder = '" << config.pString("tempSubFolder") << "'" << endl;
	cout << endl;

	// get properties for all subgroups starting with prefix
	map<string, Config*> messages = config.getGroups(); // all groups
	const string messagePrefix = "message"; // prefix for group name
	for (map<string, Config*>::iterator i = messages.begin(); i != messages.end(); ++i) {
		string groupName = i->first;
		Config* group = i->second;

		// test group name for prefix
		if (groupName.substr(0, messagePrefix.length()) == messagePrefix) {
			// display group contents
			cout << group->pString("name") << ":" << endl;
			cout << "   " << group->pString("text") << endl;
		}
	}

	return 0;
}

