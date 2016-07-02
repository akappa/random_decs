#ifndef __IOS_
#define __IOS_

#include <stdexcept>
#include <fstream>
#include <string>
#include <memory>

class ioexception : public std::runtime_error {
public:
    explicit ioexception(const char *what) : std::runtime_error(what) { }
};

void open_file(std::ifstream &file, const char *name);


void open_file(std::ofstream &file, const char *name);

std::streamoff file_length(std::ifstream &f) throw (ioexception);

template <typename T>
void read_some(std::ifstream &f, T *data, std::streamsize length) throw (ioexception)
{
    f.read(reinterpret_cast<char*>(data), sizeof(T) * length);
}

template <typename T>
void read_data(std::ifstream &f, T *data, std::streamsize length) throw (ioexception)
{
    f.read(reinterpret_cast<char*>(data), sizeof(T) * length);
    if (f.bad() || f.fail())
        throw ioexception("Failed to read file");
}

template <typename T>
void read_file(std::ifstream &f, T *data, std::streamsize length) throw (ioexception)
{
    f.seekg(0, std::ios_base::beg);
    read_data(f, data, length);
}

template <typename T>
void write_file(std::ofstream &f, T *data, std::streamsize length) throw (ioexception)
{
    f.write(reinterpret_cast<char*>(data), sizeof(T) * length);
    if (f.bad() || f.fail())
        throw ioexception("Failed to write on file");
}

template <typename T = char>
std::unique_ptr<T[]> read_file(const char *name, size_t *size, unsigned int extra = 0)
{
    std::ifstream file;
    open_file(file, name);
    std::streamoff len = file_length(file);
    std::unique_ptr<T[]> to_ret(new T[len + extra]);
    read_file(file, to_ret.get(), static_cast<std::streamsize>(len));
    *size = len;
    return std::move(to_ret);
}

#endif
