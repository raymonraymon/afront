
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



#ifndef GTB_COMMANDLINE_H
#define GTB_COMMANDLINE_H

#include <vector>
#include <gtb/common.hpp>

GTB_BEGIN_NAMESPACE

class CommandLine {

public:
	CommandLine(int argc, char **argv, const char *usage) : _argc(argc), _argv(argv), _usage(usage) { }

	void AddVariable(const char *name, bool *var, const char *help);
	void AddVariable(const char *name, int *var, const char *help);
	void AddVariable(const char *name, float *var, const char *help);
	void AddVariable(const char *name, double *var, const char *help);
	void AddVariable(const char *name, char **var, const char *help);

    typedef void (*PFUN)();
    typedef int (*PFUN2)(int argc, char** argv);
	void AddFunction(const char *name, PFUN,  const char *help);
	void AddFunction(const char *name, PFUN2, const char *help);

    void Print();
	void Parse(int start);

private:
	int _argc;
	char **_argv;
	const char *_usage;


	enum OptionType {
		CharVariable,
		BoolVariable,
		IntVariable,
		FloatVariable,
		DoubleVariable,
		Function,
		FunctionNoArgs
	};


	struct Option {
		const char *name;
		OptionType type;
		PFUN  pfun;
		PFUN2 pfun2;
		void *varp;
		const char *help;
	};

	std::vector<Option> _options;

};

#define CL_ADD_VAR(cl,name,help) cl.AddVariable(#name, &name, help)
#define CL_ADD_FUN(cl,name,help) cl.AddFunction(#name, do_##name, help)

GTB_END_NAMESPACE

#endif
