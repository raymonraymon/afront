
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
#ifndef WIN32
#include <gtb/network/network.hpp>
#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#endif // WIN32
using namespace gtb;


void work(int client_socket)
{
	assert(client_socket >= 0);
	const char message[] = "Hello, world!";
	write_message(client_socket, message, strlen(message));
	close(client_socket);
}


int main(int argc, char *argv[])
{
	if (2 != argc) {
		fprintf(stderr, "Usage: %s <port_num>\n", argv[0]);
		exit(EXIT_FAILURE);
	}

	int port = strtol(argv[1], 0, 10);
	if ((LONG_MAX == port) || (LONG_MIN == port)) {
		perror("strtol");
		exit(EXIT_FAILURE);
	}

	socket_server server;
	server.connect(port);
	for (;;) {
		int client_socket = server.accept();
		work(client_socket);
	}
	/*NOTREACHED*/
	return 0;
}
