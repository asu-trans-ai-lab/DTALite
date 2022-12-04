
float g_read_float(FILE* f)
{
    if (feof(f) == 1)
        return -1;
    /*
        read a floating point number from the current pointer of the file,
        skip all spaces
     */
    char ch, buf[32];
    int i = 0;
    int flag = 1;

    /* returns -1 if end of file is reached */
    while (true)
    {
        ch = getc(f);
        if (ch == EOF || ch == '*' || ch == '$' || ch < -1 || ch >=255 )
            return -1;

        if (isdigit(ch))
            break;

        if (ch == '-')
            flag = -1;
        else
            flag = 1;
    }

    if (ch == EOF) return -1;
    while (isdigit(ch) || ch == '.') {
        buf[i++] = ch;
        ch = fgetc(f);

    }
    buf[i] = 0;

    /* atof function converts a character string (char *) into a doubleing
    pointer equivalent, and if the string is not a floting point number,
    a zero will be return.
    */

    return (float)(atof(buf) * flag);
}

int g_ParserIntSequence(std::string str, std::vector<int>& vect)
{
    std::stringstream ss(str);
    int i;

    while (ss >> i)
    {
        vect.push_back(i);
        if (ss.peek() == ';')
            ss.ignore();
    }

    return vect.size();
}


int g_ParserStringSequence(std::string str_input, vector<string>& vect)
{

    std::istringstream ss(str_input);
    std::string token;

    while (std::getline(ss, token, ';')) {
        vect.push_back(token);
    }

    return vect.size();
}


vector<float> g_time_parser(string str)
{
    vector<float> output_global_minute;

    int string_lenghth = str.length();

    //ASSERT(string_lenghth < 100);

    const char* string_line = str.data(); //string to char*

    int char_length = strlen(string_line);

    char ch, buf_ddhhmm[32] = { 0 }, buf_SS[32] = { 0 }, buf_sss[32] = { 0 };
    char dd1, dd2, hh1, hh2, mm1, mm2, SS1, SS2, sss1, sss2, sss3;
    float ddf1, ddf2, hhf1, hhf2, mmf1, mmf2, SSf1, SSf2, sssf1, sssf2, sssf3;
    float global_minute = 0;
    float dd = 0, hh = 0, mm = 0, SS = 0, sss = 0;
    int i = 0;
    int buffer_i = 0, buffer_k = 0, buffer_j = 0;
    int num_of_colons = 0;

    //DDHHMM:SS:sss or HHMM:SS:sss

    while (i < char_length)
    {
        ch = string_line[i++];

        if (num_of_colons == 0 && ch != '_' && ch != ':') //input to buf_ddhhmm until we meet the colon
        {
            buf_ddhhmm[buffer_i++] = ch;
        }
        else if (num_of_colons == 1 && ch != ':') //start the Second "SS"
        {
            buf_SS[buffer_k++] = ch;
        }
        else if (num_of_colons == 2 && ch != ':') //start the Millisecond "sss"
        {
            buf_sss[buffer_j++] = ch;
        }

        if (ch == '_' || i == char_length) //start a new time string
        {
            if (buffer_i == 4) //"HHMM"
            {
                //HHMM, 0123
                hh1 = buf_ddhhmm[0]; //read each first
                hh2 = buf_ddhhmm[1];
                mm1 = buf_ddhhmm[2];
                mm2 = buf_ddhhmm[3];

                hhf1 = ((float)hh1 - 48); //convert a char to a float
                hhf2 = ((float)hh2 - 48);
                mmf1 = ((float)mm1 - 48);
                mmf2 = ((float)mm2 - 48);

                dd = 0;
                hh = hhf1 * 10 * 60 + hhf2 * 60;
                mm = mmf1 * 10 + mmf2;
            }
            else if (buffer_i == 6) //"DDHHMM"
            {
                //DDHHMM, 012345
                dd1 = buf_ddhhmm[0]; //read each first
                dd2 = buf_ddhhmm[1];
                hh1 = buf_ddhhmm[2];
                hh2 = buf_ddhhmm[3];
                mm1 = buf_ddhhmm[4];
                mm2 = buf_ddhhmm[5];

                ddf1 = ((float)dd1 - 48); //convert a char to a float
                ddf2 = ((float)dd2 - 48);
                hhf1 = ((float)hh1 - 48);
                hhf2 = ((float)hh2 - 48);
                mmf1 = ((float)mm1 - 48);
                mmf2 = ((float)mm2 - 48);

                dd = ddf1 * 10 * 24 * 60 + ddf2 * 24 * 60;
                hh = hhf1 * 10 * 60 + hhf2 * 60;
                mm = mmf1 * 10 + mmf2;
            }

            if (num_of_colons == 1 || num_of_colons == 2)
            {
                //SS, 01
                SS1 = buf_SS[0]; //read each first
                SS2 = buf_SS[1];

                SSf1 = ((float)SS1 - 48); //convert a char to a float
                SSf2 = ((float)SS2 - 48);

                SS = (SSf1 * 10 + SSf2) / 60;
            }

            if (num_of_colons == 2)
            {
                //sss, 012
                sss1 = buf_sss[0]; //read each first
                sss2 = buf_sss[1];
                sss3 = buf_sss[2];

                sssf1 = ((float)sss1 - 48); //convert a char to a float
                sssf2 = ((float)sss2 - 48);
                sssf3 = ((float)sss3 - 48);

                sss = (sssf1 * 100 + sssf2 * 10 + sssf3) / 1000;
            }

            global_minute = dd + hh + mm + SS + sss;

            output_global_minute.push_back(global_minute);

            //initialize the parameters
            buffer_i = 0;
            buffer_k = 0;
            buffer_j = 0;
            num_of_colons = 0;
        }

        if (ch == ':')
            num_of_colons += 1;
    }

    return output_global_minute;
}

string g_time_coding(float time_stamp)
{
    int hour = static_cast<int>(time_stamp / 60);
    int minute = static_cast<int>(time_stamp - hour * 60);
    int second = static_cast<int>((time_stamp - hour * 60 - minute) * 60 + 0.02);

    int sss = ((time_stamp - hour * 60 - minute) * 60 - second)*1000;

    //mm:ss.sss
    ostringstream strm;
    strm.fill('0');
    strm << std::setw(2) << hour << std::setw(2) << minute << ":" << std::setw(2) << second << "." << std::setw(3) << sss;

    return strm.str();
}

