/* Portions Copyright 2019 Xuesong Zhou and Peiheng Li
 *
 * If you help write or modify the code, please also list your names here.
 * The reason of having Copyright info here is to ensure all the modified version, as a whole, under the GPL
 * and further prevent a violation of the GPL.
 *
 * More about "How to use GNU licenses for your own software"
 * http://www.gnu.org/licenses/gpl-howto.html
 *
 *
 * A class that simultaneously outputs to ostream and ofstream objects
 * source: http://www.cplusplus.com/forum/general/64174/#msg347154
 *
 * Minor changes by Peiheng Li include,
 *  1. replacing typedef with using;
 *  2. removing the empty space between std::char_traits<C> and '>' in
 *     template <typename C, typename T = std::char_traits<C>> as compliers
 *     complying C++11 and higer can differentiate the ending ">>" with the istream
 *     operator ">>".
 */

#ifndef GUARD_TEESTREAM_H
#define GUARD_TEESTREAM_H

#include <iostream>

template <typename C, typename T = std::char_traits<C>>

struct basic_teebuf : public std::basic_streambuf<C, T>
{
    using streambuf_type = std::basic_streambuf<C, T>;
    using int_type = typename T::int_type;

    basic_teebuf(streambuf_type* buff_a, streambuf_type* buff_b)
        : first(buff_a), second(buff_b) {}

protected:
    virtual int_type overflow(int_type c)
    {
        const int_type eof = T::eof();
        if (T::eq_int_type(c, eof)) return T::not_eof(c);
        else
        {
            const C ch = T::to_char_type(c);
            if (T::eq_int_type(first->sputc(ch), eof) ||
                T::eq_int_type(second->sputc(ch), eof))
                return eof;
            else return c;
        }
    }

    virtual int sync()
    {
        return !first->pubsync() && !second->pubsync() ? 0 : -1;
    }

private:
    streambuf_type* first;
    streambuf_type* second;
};

template <typename C, typename T = std::char_traits<C>>
struct basic_teestream : public std::basic_ostream<C, T>
{
    basic_teestream() : debug_level{ 0 }, log_sig{ 0 }, log_odme{ 0 }, log_path{ 0 } , log_ue{ 0 }
    {
    
    }
    // add more controls of the debug logs
    int debug_level;
    int log_sig;
    int log_odme;
    int log_path;
    int log_dta;
    int log_ue;

    using stream_type = std::basic_ostream<C, T>;
    using streambuff_type = basic_teebuf<C, T> ;

    basic_teestream(stream_type& first, stream_type& second)
        : stream_type(&stmbuf), stmbuf(first.rdbuf(), second.rdbuf()) {}

    basic_teestream(streambuff_type* first, streambuff_type* second)
        : stream_type(&stmbuf), stmbuf(first, second) {}

    ~basic_teestream() { stmbuf.pubsync(); }

private:
    streambuff_type stmbuf;
};

using teebuf = basic_teebuf<char>;
using teestream = basic_teestream<char>;

#endif