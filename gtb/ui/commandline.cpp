
/*

Copyright 2007 University of Utah


This file is part of Afront.

Afront is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published
by the Free Software Foundation; either version 2 of the License,
or (at your option) any later version.

Afront is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA

*/



#include <gtb/gtb.hpp>

#include <vector>
#include <iostream>
#include <string>
#include "commandline.h"

#ifndef WIN32
#define stricmp strcasecmp
#endif

GTB_BEGIN_NAMESPACE

void CommandLine::AddVariable(const char *name, char **var, const char *help) {
	struct Option opt;
	opt.name = name;
	opt.type = CharVariable;
	opt.varp = var;
	opt.help = help;
	_options.push_back(opt);
}

void CommandLine::AddVariable(const char *name, bool *var, const char *help) {
	struct Option opt;
	opt.name = name;
	opt.type = BoolVariable;
	opt.varp = var;
	opt.help = help;
	_options.push_back(opt);
}

void CommandLine::AddVariable(const char *name, int *var, const char *help) {
	struct Option opt;
	opt.name = name;
	opt.type = IntVariable;
	opt.varp = var;
	opt.help = help;
	_options.push_back(opt);
}

void CommandLine::AddVariable(const char *name, float *var, const char *help) {
	struct Option opt;
	opt.name = name;
	opt.type = FloatVariable;
	opt.varp = var;
	opt.help = help;
	_options.push_back(opt);
}

void CommandLine::AddVariable(const char *name, double *var, const char *help) {
	struct Option opt;
	opt.name = name;
	opt.type = DoubleVariable;
	opt.varp = var;
	opt.help = help;
	_options.push_back(opt);
}


void CommandLine::AddFunction(const char *name, PFUN fun, const char *help) {
	struct Option opt;
	opt.name = name;
	opt.type = FunctionNoArgs;
	opt.pfun = fun;
	opt.help = help;
	_options.push_back(opt);
}


void CommandLine::AddFunction(const char *name, PFUN2 fun, const char *help) {
	struct Option opt;
	opt.name = name;
	opt.type = Function;
	opt.pfun2 = fun;
	opt.help = help;
	_options.push_back(opt);
}


void CommandLine::Print() {

	std::cerr<<"usage: "<<_usage<<std::endl;
	for (unsigned i=0; i<_options.size(); i++) {

		std::string str = "-";
		str += _options[i].name;

		const char *colon = strchr(_options[i].help, ':');
		if (colon) {
			std::string args = _options[i].help;
			args.erase(colon-_options[i].help, strlen(colon));
			str += " ";
			str += args;
			while (str.size()<30) str += " ";
			str += colon;
		} else {
			while (str.size()<30) str += " ";
			str += _options[i].help;
		}

		std::cerr<<str;


		switch (_options[i].type) {
		    
		case BoolVariable:
		    std::cerr<<"     =["<<(*(bool*)_options[i].varp ? "t" : "f")<<"]";
		    break;
		case CharVariable:
		    if ((*(char**)_options[i].varp) == NULL)
			std::cerr<<"     =[NULL]";
		    else
			std::cerr<<"     =["<<(*(char**)_options[i].varp)<<"]";
		    break;
		case IntVariable:
		    std::cerr<<"     =["<<(*(int*)_options[i].varp)<<"]";
		    break;
		case FloatVariable:
		    std::cerr<<"     =["<<(*(float*)_options[i].varp)<<"]";
		    break;
		case DoubleVariable:
		    std::cerr<<"     =["<<(*(double*)_options[i].varp)<<"]";
		    break;
		case Function:
		    break;
		case FunctionNoArgs:
		    break;
		default:
		    exit(-1);
		}

		 
		std::cerr<<std::endl;
	}
}

void CommandLine::Parse(int start) {
	int cur=start;
	while (cur<_argc) {

		if (_argv[cur][0] != '-' || !stricmp(_argv[cur], "-?")) {
			Print(); exit(-1);
		}

		bool found=false;
		for (unsigned i=0; i<_options.size(); i++) {

			if (!stricmp(_options[i].name, _argv[cur]+1)) {
				found=true;
				switch (_options[i].type) {

					case BoolVariable:
						if (cur+1 >= _argc) { Print(); exit(-1); }
						if (!stricmp(_argv[cur+1], "true") || !stricmp(_argv[cur+1], "t"))
							*(bool*)_options[i].varp = true;
						else if (!stricmp(_argv[cur+1], "false") || !stricmp(_argv[cur+1], "f"))
							*(bool*)_options[i].varp = false;
						else
							*(bool*)_options[i].varp = (atoi(_argv[cur+1])!=0) ? true : false;
						cur+=2;
						break;
					case CharVariable:
						if (cur+1 >= _argc) { Print(); exit(-1); }
						*(char**)_options[i].varp = _argv[cur+1];
						cur+=2;
						break;
					case IntVariable:
						if (cur+1 >= _argc) { Print(); exit(-1); }
						*(int*)_options[i].varp = atoi(_argv[cur+1]);
						cur+=2;
						break;
					case FloatVariable:
						if (cur+1 >= _argc) { Print(); exit(-1); }
						*(float*)_options[i].varp = (float)atof(_argv[cur+1]);
						cur+=2;
						break;
					case DoubleVariable:
						if (cur+1 >= _argc) { Print(); exit(-1); }
						*(double*)_options[i].varp = atof(_argv[cur+1]);
						cur+=2;
						break;
					case Function:
						cur += _options[i].pfun2(_argc-cur, _argv+cur);
						break;
					case FunctionNoArgs:
						_options[i].pfun();
						cur++;
						break;
					default:
						exit(-1);

				}

				break;
			}

			if (i == _options.size()) { Print(); exit(-1); }
		}

		if (!found) { Print(); exit(-1); }

	}

}

GTB_END_NAMESPACE







