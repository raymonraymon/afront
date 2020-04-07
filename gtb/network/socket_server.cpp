
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
#include <gtb/network/socket_server.hpp>
#include <gtb/network/netutil.hpp>
//#include <config.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <arpa/inet.h>
#include <unistd.h>
#endif // WIN32


GTB_BEGIN_NAMESPACE


socket_server::socket_server()
{
}


bool socket_server::connect(int port)
{
	// create the server socket
	_server_socket = socket(PF_INET, SOCK_STREAM, IPPROTO_TCP);
	if (-1 == _server_socket) {
		perror("socket");
		//exit(EXIT_FAILURE);
		close(_server_socket);
		return false;
	}
	fprintf(stderr, "server socket: %d\n", _server_socket);

	// reuse the address quickly if the server dies
	int on = 0;
	setsockopt(_server_socket, SOL_SOCKET, SO_REUSEADDR,
		   (const char *) &on, sizeof(on));

	// linger on sending data once the connection closes
	struct linger my_linger;
  	my_linger.l_onoff = 1;
  	my_linger.l_linger = 30;
	setsockopt(_server_socket, SOL_SOCKET, SO_LINGER,
		   (const char *) &my_linger, sizeof(my_linger));

	// get hostname
	char host_name[80];
	if (-1 == get_host_name(host_name, sizeof(host_name))) {
		perror("get_host_name");
		//exit(EXIT_FAILURE);
		close(_server_socket);
		return false;
	}
	fprintf(stderr, "server host name: %s\n", host_name);

	struct hostent *host_ptr = 0;
	host_ptr = gethostbyname(host_name);
	if (0 == host_ptr) {
		perror("gethostbyname");
		//exit(EXIT_FAILURE);
		close(_server_socket);
		return false;
	}

	// create and bind to ip address:port
	struct sockaddr_in server_name;
	memset(&server_name, 0, sizeof(server_name));
	memcpy(&server_name.sin_addr, host_ptr->h_addr, host_ptr->h_length);

	server_name.sin_addr.s_addr = htonl(INADDR_ANY);
	server_name.sin_family = AF_INET;
	server_name.sin_port = htons(port);

	if (-1 == bind(_server_socket, (struct sockaddr *) &server_name,
		       sizeof(server_name))) {
		perror("bind");
		//exit(EXIT_FAILURE);
		close(_server_socket);
		return false;
	}

	// listen for incoming connection requests
	if (-1 == listen(_server_socket, 5)) {
		perror("listen");
		//exit(EXIT_FAILURE);
		close(_server_socket);
		return false;
	}
	return true;
}


socket_server::~socket_server()
{
	close(_server_socket);
}


int socket_server::accept()
{
	struct sockaddr_in client_name;
	socklen_t client_len = sizeof(client_name);

	// accept the connection
	memset(&client_name, 0, sizeof(client_name));
	int client_socket = ::accept(_server_socket,
				     (struct sockaddr *) &client_name,
				     &client_len);
	if (-1 == client_socket) {
		perror("accept");
		exit(EXIT_FAILURE);
	}

	// find out who the client peer is
	if (-1 == getpeername(client_socket,
			      (struct sockaddr *) &client_name,
			      &client_len)) {
		perror("getpeername");
		//exit(EXIT_FAILURE);
	} else {
		printf("Connection request from %s\n",
		       inet_ntoa(client_name.sin_addr));
	}
	return client_socket;
}


GTB_END_NAMESPACE
