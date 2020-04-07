
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
#include <gtb/network/netutil.hpp>
#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/utsname.h>
#include <netinet/in.h>		// for hton* and ntoh*
#endif // WIN32


int get_host_name(char *buffer, int len)
{
#ifdef WIN32
	DWORD nSize = len;
	BOOL rc = GetComputerName(buffer, &nSize);
	if (rc) return 0;
	else return -1;
#else // WIN32
	struct utsname sys_name;
	int status = uname(&sys_name);
	if (-1 != status) {
		strncpy(buffer, sys_name.nodename, len);
	}
	return status;
#endif // WIN32
}


void read_message(int sock_descr, void *msg_ptr, int bytes_to_read)
{
	char *msg = (char *) msg_ptr;
	for (int bytes_read = 0, n = 0;
	     bytes_read < bytes_to_read;
	     bytes_read += n) {
		n = read(sock_descr,
			 msg + bytes_read,
			 (bytes_to_read - bytes_read));
		if (n < 0) {
			perror("read");
			if (errno != EINTR) {
				exit(EXIT_FAILURE);
			}
		}
	}
}


void write_message(int sock_descr, const void *msg_ptr, int bytes_to_write)
{
	const char *msg = (const char *) msg_ptr;
	for (int bytes_written = 0, n = 0;
	     bytes_written < bytes_to_write;
	     bytes_written += n) {
		n = write(sock_descr,
			  msg + bytes_written,
			  (bytes_to_write - bytes_written));
		if (n < 0) {
			perror("write");
			if (errno != EINTR) {
				exit(EXIT_FAILURE);
			}
		}
	}
}


// reads

void net_read_bool(bool &b, int fd)
{
	read_message(fd, &b, sizeof(b));
}


void net_read_char(char &c, int fd)
{
	read_message(fd, &c, sizeof(c));
}


void net_read_uchar(unsigned char &c, int fd)
{
	read_message(fd, &c, sizeof(c));
}


void net_read_int(int &x, int fd)
{
	read_message(fd, &x, sizeof(x));
	x = ntohl(x);
}


void net_read_uint(unsigned &x, int fd)
{
	read_message(fd, &x, sizeof(x));
	x = ntohl(x);
}


void net_read_float(float &x, int fd)
{
	unsigned i;
	net_read_uint(i, fd);
	x = *(float *) &i;
}


void net_read_str(char *buffer, unsigned size, int fd)
{
	unsigned n;
	net_read_uint(n, fd);
	assert(n < size);
	(void) size;
	read_message(fd, buffer, n);
	buffer[n] = '\0';
}


// writes

void net_write_bool(bool b, int fd)
{
	write_message(fd, &b, sizeof(b));
}


void net_write_char(char c, int fd)
{
	write_message(fd, &c, sizeof(c));
}


void net_write_uchar(unsigned char c, int fd)
{
	write_message(fd, &c, sizeof(c));
}


void net_write_int(int x, int fd)
{
	x = htonl(x);
	write_message(fd, &x, sizeof(x));
}


void net_write_uint(unsigned x, int fd)
{
	x = htonl(x);
	write_message(fd, &x, sizeof(x));
}


void net_write_float(float x, int fd)
{
	unsigned i = *(unsigned *) &x;
	net_write_uint(i, fd);
}


void net_write_str(const char *s, int fd)
{
	unsigned n = strlen(s);
	net_write_uint(n, fd);
	write_message(fd, s, n);
}
