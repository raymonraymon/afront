
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


#ifndef GTB_SOCKET_CLIENT_INCLUDED
#define GTB_SOCKET_CLIENT_INCLUDED

#include <gtb/common.hpp>

GTB_BEGIN_NAMESPACE

class socket_client {
public:
	socket_client();

	~socket_client();

	bool connect(const char *remote_host, int remote_port);

	int socket();

protected:
	int _socket;
};

GTB_END_NAMESPACE

#endif // GTB_SOCKET_CLIENT_INCLUDED
