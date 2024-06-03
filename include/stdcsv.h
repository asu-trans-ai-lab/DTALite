/**
 * @file stdcsv.h, part of the project MIOCSV under Apache License 2.0
 * @author jdlph (jdlph@hotmail.com)
 * @brief A suite of lightening-fast CSV parsers and writer built upon the C++ std facilities
 *
 * @copyright Copyright (c) 2022 - 2023 Peiheng Li, Ph.D.
 */

#ifndef GUARD_STDCSV_H
#define GUARD_STDCSV_H

#define O3N_TIME_BOUND

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace miocsv
{
using size_type = unsigned long;
using FieldNames = std::map<std::string, size_type>;

/**
 * @brief a helper class to define a string range by [head, tail]
 *
 * it is to mimic std::string_view but with capability in changing head and tail.
 * we use it to avoid potential dynamic memory allocation (and deallocation) and
 * copy in string concatenation.
 */
template <typename InputIter>
class StringRange {
public:
    StringRange() = delete;

    explicit StringRange(InputIter it) : head {it}, tail {it}
    {
    }

    StringRange(const StringRange&) = delete;
    StringRange& operator=(const StringRange&) = delete;

    StringRange(StringRange&&) = delete;
    StringRange& operator=(StringRange&&) = delete;

    ~StringRange() = default;

    char back() const
    {
        return *(tail - 1);
    }

    bool empty() const
    {
        return tail == head;
    }

    void extend(InputIter it)
    {
        // check it > tail?
        tail = it;
    }

    void reset(InputIter it)
    {
        head = tail = it;
    }

    std::string to_string() const
    {
        return std::string{head, tail};
    }

    // in case that the last char is '\r'
    std::string to_string_cr() const
    {
        return std::string{head, tail - 1};
    }

private:
    InputIter head;
    InputIter tail;
};

/**
 * @brief a custom exception which arises when a nonexistent field is retrieved
 * via Row::[] (i.e., no such a field in the header)
 */
struct NoRecord : public std::runtime_error {
    NoRecord() = delete;

    explicit NoRecord(const std::string& s)
        : std::runtime_error{"Row::operator[]: " + s + " is nonexistent"}
    {
    }

    explicit NoRecord(size_type s)
        : std::runtime_error{"Row::operator[]: " + std::to_string(s) + " is out of range"}
    {
    }
};

class Row {
    friend void attach_fieldnames(Row&, const FieldNames*, size_type);
    friend std::ostream& operator<<(std::ostream&, const Row&);

public:
    using Records = std::vector<std::string>;
    using iterator = Records::iterator;
    using const_iterator = Records::const_iterator;

    Row() = default;

    // pure string input
    Row(std::initializer_list<std::string> args)
    {
        for (auto& s: args)
            records.emplace_back(s);
    }

    template<typename T, typename... Args>
    Row(const T& t, const Args&... args)
    {
        convert_to_string(t, args...);
    }

    Row(Records&& r)
    {
        records = r;
    }

    Row(const Row&) = default;
    Row& operator=(const Row&) = delete;

    Row(Row&&) = default;
    Row& operator=(Row&&) = default;

    ~Row() = default;

    std::string& operator[](size_type i)
    {
        if (i < 0 || i >= records.size())
            throw NoRecord{i};

        return records[i];
    }

    const std::string& operator[](size_type i) const
    {
        if (i < 0 || i >= records.size())
            throw NoRecord{i};

        return records[i];
    }

    // use std::runtime_error?
    std::string& operator[](const std::string& s)
    {
        try
        {
            size_type i = fns->at(s);
            // more fieldnames than fields will be taken care by operator[]
            return records[i];
        }
        catch (const std::out_of_range)
        {
            throw NoRecord{s};
        }
    }

    /**
     * @brief retrieve a field (record) using fieldname (header)
     *
     * in case of data inconsistency, the following two cases are considered.
     *
     * 1. if there are more fields than fieldnames, users can only retrieve those
     * extra fields (with no corresponding fieldnames) using indices via overloaded
     * [].
     * 2. if there are more fieldnames than fields, an exception NoRecord will be
     * thrown at run time.
     *
     * Note that our way handling these two exceptions are different with Python
     * csv.DictReader, where, extra fields will be concatenated as a list
     * (of strings) with None assigned as a new fieldname for case 1; None will
     * be output for case 2.
     *
     * @param s
     * @return const std::string&
     */
    const std::string& operator[](const std::string& s) const
    {
        try
        {
            size_type i = fns->at(s);
            // more fieldnames than fields will be taken care by operator[]
            return records[i];
        }
        catch (const std::out_of_range)
        {
            throw NoRecord{s};
        }
    }

    std::string& back()
    {
        return records.back();
    }

    const std::string& back() const
    {
        return records.back();
    }

    iterator begin()
    {
        return records.begin();
    }

    const_iterator begin() const
    {
        return records.begin();
    }

    iterator end()
    {
        return records.end();
    }

    const_iterator end() const
    {
        return records.end();
    }

    size_type size() const
    {
        return records.size();
    }

    bool empty() const
    {
        return records.size() == 0;
    }

    // move non-const string
    void append(std::string& s)
    {
        records.emplace_back(s);
    }

    // move rvalue string
    void append(std::string&& s)
    {
        records.push_back(s);
    }

    // copy const string with no move
    void append(const std::string& s)
    {
        records.push_back(s);
    }

    void append(const std::string_view sv)
    {
        records.push_back(std::string{sv});
    }

private:
    Records records;
    // reserved for DictReader
    const FieldNames* fns = nullptr;

    template<typename T>
    void convert_to_string(const T& t)
    {
        std::ostringstream os;
        os << t;
        // os.str() will be moved into records as os.str() is a rvalue reference
        records.push_back(os.str());
    }

    template<typename T, typename... Args>
    void convert_to_string(const T& t, const Args&... args)
    {
#if __cplusplus >= 201703L
        std::ostringstream os;
        os << t;
        records.push_back(os.str());
        // it requires C++17
        if constexpr(sizeof...(args) > 0)
            convert_to_string(args...);
#else
        convert_to_string(t);
        convert_to_string(args...);
#endif
    }
};

// not a pure abstract class
class BaseReader {
public:
    BaseReader() : quote {'"'}, row_num {0}
    {
    }

    BaseReader(const BaseReader&) = delete;
    BaseReader& operator=(const BaseReader&) = delete;

    BaseReader(BaseReader&&) = default;
    BaseReader& operator=(BaseReader&&) = delete;

    virtual ~BaseReader()
    {
    }

    // forward declaration
    class ReaderIterator;

    ReaderIterator begin();
    ReaderIterator end();

    size_type get_row_num() const
    {
        return row_num;
    }

protected:
    const char quote;
    size_type row_num;

    Row row;

    const char CR = '\r';
    const char LF = '\n';

    struct InvalidRow : public std::runtime_error {
        // do we need it?
        InvalidRow() = delete;

        InvalidRow(size_type row_num, const std::string& str)
            : std::runtime_error{"CAUTION: Invalid Row at line "
                                 + std::to_string(++row_num)
                                 + "! Value is not allowed after quoted field: "
                                 + str}
        {
        }
    };

    class IterationEnd {

    };

    virtual void iterate() = 0;
};

class BaseReader::ReaderIterator {
public:
    ReaderIterator() = delete;

    ReaderIterator(BaseReader* r_) : r {r_}
    {
        if (!r)
            return;

        try
        {
            r->iterate();
        }
        catch (IterationEnd)
        {
            r = nullptr;
        }
    }

    ReaderIterator(const ReaderIterator&) = default;
    ReaderIterator& operator=(const ReaderIterator&) = delete;

    ReaderIterator(ReaderIterator&&) = default;
    ReaderIterator& operator=(ReaderIterator&&) = delete;

    ~ReaderIterator() = default;

    ReaderIterator operator++()
    {
        try
        {
            r->iterate();
        }
        catch (IterationEnd)
        {
            r = nullptr;
        }

        return *this;
    }

    bool operator==(const ReaderIterator& it) const
    {
        return r == it.r;
    }

    bool operator!=(const ReaderIterator& it) const
    {
        return r != it.r;
    }

    const Row& operator*() const
    {
        return r->row;
    }

private:
    BaseReader* r;
};

// dreaded diamond
class BaseDictReader : public virtual BaseReader {
public:
    BaseDictReader() = default;

    void setup_headers(const Row& r)
    {
        for (size_type i = 0, sz = r.size(); i != sz; ++i)
        {
            const auto& s = r[i];
            fns[s] = i;
        }

        if (fns.empty() && row_num == 0)
        {
            try
            {
                iterate();
                setup_headers(row);
            }
            catch (IterationEnd)
            {
                return;
            }
        }
    }

    const FieldNames& get_fieldnames() const
    {
        return fns;
    }

protected:
    FieldNames fns;
};

// dreaded diamond
class Reader : public virtual BaseReader {
public:
    Reader() = delete;

    Reader(const std::string& ist_, const char delim_ = ',')
        : BaseReader{}, ist {ist_}, delim {delim_}
    {
        if (!ist)
        {
            std::cerr << "invalid input! no " << ist_ << '\n';
            std::terminate();
        }
#ifdef O3N_TIME_BOUND
        it = ist;
#endif
    }

    Reader(std::string&& ist_, const char delim_ = ',')
        : BaseReader{}, ist {ist_}, delim {delim_}
    {
        if (!ist)
        {
            std::cerr << "invalid input! no " << ist_ << '\n';
            std::terminate();
        }
#ifdef O3N_TIME_BOUND
        it = ist;
#endif
    }

protected:
    std::ifstream ist;
    const char delim;

    void iterate() override
    {
#ifdef O3N_TIME_BOUND
        if (it == it_end)
            throw IterationEnd{};
#else
        std::string s;
        if (!std::getline(ist, s))
            throw IterationEnd{};
#endif

        try
        {
#ifdef O3N_TIME_BOUND
            row = split3();
#else
            row = split2(s);
#endif
        }
        catch (const InvalidRow& e)
        {
            std::cerr << e.what() << '\n';
            // this is necessary. otherwise, DictReader::iterate() will be called
            // if DictReader is being used.
            Reader::iterate();
        }

        ++row_num;
    }

private:
#ifdef O3N_TIME_BOUND
    std::istreambuf_iterator<char> it;
    std::istreambuf_iterator<char> it_end;
#endif

    // for benchmark only
    Row split(const std::string& s) const;

    /**
     * @brief parse string
     *
     * @param C string container, which could be std::string or std::string_view
     * @return Row
     *
     * @note with split2(), the overall time complexity is O(5N) which includes two
     * linear searches and three copy processes.
     */
    template<typename C>
    Row split2(const C& c) const;

    /**
     * @brief parse string
     *
     * @return Row
     *
     * @note with split3(), the overall time complexity is O(3N) which includes one
     * linear search and two copy processes.
     */
    Row split3();
};

class DictReader : public Reader, public BaseDictReader {
public:
    DictReader() = delete;

    DictReader(const std::string& ist_, const Row& fieldnames_ = {}, const char delim_ = ',')
        : Reader{ist_}, BaseDictReader{}
    {
        setup_headers(fieldnames_);
    }

    DictReader(std::string&& ist_, const Row& fieldnames_ = {}, const char delim_ = ',')
        : Reader{ist_}, BaseDictReader{}
    {
        setup_headers(fieldnames_);
    }

private:
    void iterate() override
    {
        Reader::iterate();
        // do not take blank lines in consistent with Python csv.DictReader
        while (row.empty())
            Reader::iterate();

        if (row_num > 1)
            attach_fieldnames(row, &fns, row_num);
    }
};

class Writer {
public:
    Writer() = delete;

    Writer(const std::string& ost_, const char delim_ = ',')
        : ost {ost_}, delim {delim_}
    {
        if (!ost)
        {
            std::cerr << "invalid input! no " << ost_ << '\n';
            std::terminate();
        }
    }

    Writer(std::string&& ost_, const char delim_ = ',')
        : ost {ost_}, delim {delim_}
    {
        if (!ost)
        {
            std::cerr << "invalid input! no " << ost_ << '\n';
            std::terminate();
        }
    }

    Writer(const Writer&) = delete;
    Writer& operator=(const Writer) = delete;

    Writer(Writer&&) = default;
    Writer& operator=(Writer&&) = delete;

    ~Writer() = default;

    /**
     * @brief append context into the file
     *
     * @note need to make sure that the context has NO delimiter
     *
     * @param t a const reference of typename T
     * @param sep separator, which could be any valid char. its default is ',' as
     * for csv files.
     */
    template<typename T>
    void append(const T& t, const char sep = ',')
    {
        ost << t << sep;
    }

    template<typename T>
    void append(T&& t, const char sep = ',')
    {
        ost << t << sep;
    }

    /**
     * @brief append context into the file
     *
     * @note need to make sure that the context has NO delimiter
     *
     * @param t a const reference of typename T
     * @param str any valid string.
     */
    template<typename T>
    void append(const T& t, const std::string& str)
    {
        ost << t << str;
    }

    template<typename T>
    void append(T&& t, const std::string& str)
    {
        ost << t << str;
    }

    /**
     * @brief write a row of records into the file
     *
     * @note if no records contain the delimiter, then a simple implementation via the
     * overloaded operator<< for Row would work fine, i.e., os << r << '\n'. It
     * has been implemented as write_row_raw()
     *
     * @param r a const reference of miocsv::Row
     */
    void write_row(const Row& r)
    {
        for (auto it = r.begin(), it_end = r.end() - 1; it != it_end; ++it)
        {
            // check if the current record contains the delimiter or not
            if (std::find(it->begin(), it->end(), delim) != it->end())
                ost << '"' << *it << '"' << delim;
            else
                ost << *it << delim;
        }

        // the last record
        if (std::find(r.back().begin(), r.back().end(), delim) != r.back().end())
            ost << '"' << r.back() << '"';
        else
            ost << r.back();

        ost << '\n';
    }

    /**
     * @brief write a row of records into the file
     *
     * @note if no records contain the delimiter, then a simple implementation via the
     * overloaded operator<< for Row would work fine, i.e., os << r << '\n'. It
     * has been implemented as write_row_raw()
     *
     * @param r a rvalue reference of miocsv::Row
     */
    void write_row(Row&& r)
    {
        for (auto it = r.begin(), it_end = r.end() - 1; it != it_end; ++it)
        {
            // check if the current record contains the delimiter or not
            if (std::find(it->begin(), it->end(), delim) != it->end())
                ost << '"' << *it << '"' << delim;
            else
                ost << *it << delim;
        }

        // the last record
        if (std::find(r.back().begin(), r.back().end(), delim) != r.back().end())
            ost << '"' << r.back() << '"';
        else
            ost << r.back();

        ost << '\n';
    }

    /**
     * @brief write a row of records having no delimiters into the file
     *
     * @note users need to make sure that each record has NO delimiter. it does
     * not have handling on delimiter inside a cell as write_row(). It simply
     * writes as is.
     *
     * @param t a const reference of typename T
     */
    template<typename T>
    void write_row_raw(const T& t)
    {
        ost << t << '\n';
    }

    /**
     * @brief write a row of records having no delimiters into the file
     *
     * @note users need to make sure that each record has NO delimiter. it does
     * not have handling on delimiter inside a cell as write_row(). It simply
     * writes as is.
     *
     * @param t a const reference of typename T
     * @param args const reference of variadic template T
     */
    template<typename T, typename... Args>
    void write_row_raw(const T& t, const Args&... args)
    {
        ost << t;
        append_cell(args...);
        ost << '\n';
    }

private:
    /**
     * @brief helper function to facilitate write_row_raw()
     *
     * @param t a const reference of typename T
     */
    template<typename T>
    void append_cell(const T& t)
    {
        ost << delim << t;
    }

    /**
     * @brief helper function to facilitate write_row_raw()
     *
     * @param t a const reference of typename T
     * @param args const reference of variadic template T
     */
    template<typename T, typename... Args>
    void append_cell(const T& t, const Args&... args)
    {
        append_cell(t);
        append_cell(args...);
    }

private:
    std::ofstream ost;
    const char delim;
};

// implementations
void attach_fieldnames(Row& r, const FieldNames* fns, size_type row_num)
{
    r.fns = fns;
    if (r.fns->size() != r.size())
    {
        std::cout << "CAUTION: Data Inconsistency at line " << row_num
                  << ": " << r.fns->size() << " fieldnames vs. "
                  << r.size() << " fields\n";
    }
}

std::ostream& operator<<(std::ostream& os, const Row& r)
{
    if (r.empty())
        return os;

    for (size_type i = 0, sz = r.size(); i != sz - 1; ++i)
        os << r.records[i] << ',';
    // last one
    os << r.back();

    return os;
}

inline BaseReader::ReaderIterator BaseReader::begin()
{
    // just in case users retrieve it after iteration starts
    if (this->row_num > 1)
        return nullptr;

    return ReaderIterator{this};
}

inline BaseReader::ReaderIterator BaseReader::end()
{
    return nullptr;
}

Row Reader::split(const std::string& s) const
{
    if (s.empty())
        return Row{};

    Row r;
    std::string s1;

    auto quoted = false;
    auto b = s.begin();
    for (auto i = s.begin(), e = s.end(); i != e; ++i)
    {
        if (*i == quote)
        {
            quoted ^= true;
            if (!quoted)
            {
                s1 += std::string(b, i + 1);
                b = i + 1;
                if (*b != quote && *b != delim && b != e)
                    throw InvalidRow{row_num, s1};
            }
        }
        else if (*i == delim && !quoted)
        {
            if (!s1.empty())
            {
                r.append(s1);
                s1.clear();
            }
            else
                r.append(std::string{b, i});

            b = i + 1;
        }
    }

    // last one
    if (!s1.empty())
    {
        if (s1.back() != CR)
            r.append(s1);
        else
            r.append(std::string{s1.begin(), --s1.end()});
    }
    else
    {
        if (s.back() != CR)
            r.append(std::string{b, s.end()});
        else
            r.append(std::string{b, --s.end()});
    }

    // use move constructor to avoid copy
    return r;
}

template<typename C>
Row Reader::split2(const C& c) const
{
    Row r;
    auto quoted = false;
    StringRange<typename C::const_iterator> sr{c.begin()};

    for (auto i = c.begin(), e = c.end();;)
    {
        if (*i == quote)
        {
            sr.extend(++i);
            quoted ^= true;
            if (!quoted && *i != quote && *i != delim && i != e)
            {
                throw InvalidRow{row_num, sr.to_string()};
            }
        }
        else if (*i == delim && !quoted)
        {
            r.append(sr.to_string());
            sr.reset(++i);
        }
        else if (i == e)
        {
            // last one
            if (sr.back() != CR)
                r.append(sr.to_string());
            else
                r.append(sr.to_string_cr());

            return r;
        }
        else
            sr.extend(++i);
    }
}

#ifdef O3N_TIME_BOUND
Row Reader::split3()
{
    Row r;
    std::string s;
    auto quoted = false;

    // caution: the last line might be null terminated rather than '\n'
    while (true)
    {
        if (*it == quote)
        {
            s.push_back(*it++);
            quoted ^= true;
            if (!quoted && *it != quote && *it != delim && *it != LF)
            {
                ++it = std::find(it, it_end, LF);
                throw Reader::InvalidRow{row_num, s};
            }
        }
        else if (*it == delim && !quoted)
        {
            r.append(s);
            s.clear();
            ++it;
        }
        else if (*it == LF)
        {
            // last one
            if (s.back() != CR)
                r.append(s);
            else
                r.append(std::string{s.begin(), --s.end()});

            ++it;
            return r;
        }
        else if (it == it_end)
            return r;
        else
            s.push_back(*it++);
    }
}
#endif

/**
 * @brief a utility function to split string according to the given delimiter
 *
 * @param c string container, which could be std::string or std::string_view
 * @param delim a single character delimiter, such as ',', ';', and so on
 */
template<typename C>
Row split(const C& c, const char delim = ',')
{
    static constexpr char quote = '"';

    Row r;
    auto quoted = false;
    StringRange<typename C::const_iterator> sr{c.begin()};

    for (auto i = c.begin(), e = c.end();;)
    {
        if (*i == quote)
        {
            sr.extend(++i);
            quoted ^= true;
        }
        else if (*i == delim && !quoted)
        {
            r.append(sr.to_string());
            sr.reset(++i);
        }
        else if (i == e)
        {
            // last one
            r.append(sr.to_string());
            return r;
        }
        else
            sr.extend(++i);
    }
}

} // namespace miocsv

std::ostream& operator<<(std::ostream& os, const miocsv::FieldNames& fns)
{
    if (fns.empty())
        return os;

    using name_pos = std::pair<std::string, miocsv::size_type>;

    std::vector<name_pos> vec {fns.begin(), fns.end()};
    // sort vec to restore the insertion order
    std::sort(vec.begin(), vec.end(), [](const name_pos& left, const name_pos& right) {
        return left.second < right.second;
    });

    for (miocsv::size_type i = 0, sz = vec.size() - 1; i != sz; ++i)
        os << vec[i].first << ',';
    // last one
    os << vec.back().first;

    return os;
}

#endif