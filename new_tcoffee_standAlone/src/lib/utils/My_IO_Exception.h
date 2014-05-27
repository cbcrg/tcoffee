
#ifndef MY_IO_EXCEPTION_H_
#define MY_IO_EXCEPTION_H_


class My_IO_Exception : public std::exception
{
public:
private:
	const char *_message;

public:
	virtual const char* what() const throw()
	{
		return _message;
	}

	My_IO_Exception(const char *message):_message(message)
	{}

	virtual ~My_IO_Exception() throw()
	{}
};




#endif /* MY_IO_EXCEPTION_H_ */
