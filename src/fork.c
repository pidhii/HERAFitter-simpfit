#include <unistd.h>
#include <sys/wait.h>
#include <string.h>
#include <stdio.h>
#include <signal.h>

extern void hf_stop_();


int cfork_(int* fd)
{
	int pipefd[2];
	if(pipe(pipefd) == -1) {
		perror("pipe");
		hf_stop_();
	}
	
	int pid = fork();
	if(pid < 0)
	{
		perror("fork");
		hf_stop_();
	}
	else if(pid == 0)
	{
		close(pipefd[0]);
		*fd = pipefd[1];
	}
	else
	{
		close(pipefd[1]);
		*fd = pipefd[0];
	}

	return pid;
}

void cwaitpid_(int* pid)
{
	if(waitpid(*(pid_t*)pid, NULL, 0) == -1)
		perror("waitpid");
}

void cwrite_double_(int* fd, double* val)
{
	if(write(*fd, val, sizeof(double)) < 0) {
		perror("write");
		hf_stop_();
	}
}

void cread_double_(int* fd, double* val)
{
	if(read(*fd, val, sizeof(double)) < 0) {
		perror("read");
		hf_stop_();
	}
}

void cexit_(int* code)
{
	_exit(*code);
}

void ckill_(int* pid)
{
	kill(*pid, SIGKILL);
}

void cclose_(int* fd)
{
	if(close(*fd) == -1)
		perror("close");
}

