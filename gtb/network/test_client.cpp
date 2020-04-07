
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
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#endif // WIN32
using namespace gtb;


void work(int client_socket)
{
	char buffer[256] = "";
	int status;
	while (0 < (status = read(client_socket,
				  buffer, sizeof(buffer) - 1))) {
		buffer[status] = '\0';
		printf("socket client: received %d bytes - \"%s\"\n",
		       status, buffer);
	}
	if (-1 == status) {
		perror("read");
	}
}


int main(int argc, char *argv[])
{
	if (3 != argc) {
		fprintf(stderr, "Usage: %s <remote_host> <remote_port>\n",
			argv[0]);
		exit(EXIT_FAILURE);
	}

	char *remote_host = argv[1];
	int remote_port = atoi(argv[2]);

	socket_client client;
	client.connect(remote_host, remote_port);
	work(client.socket());

	return 0;
}
