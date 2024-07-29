#include "similarity.h"

/*
 Simplified from:
 https://github.com/wernsey/miscsrc/blob/master/simil.h
*/


int rsimil(const char* a, int alen, const char* b, int blen)
{    
    int i, j, k, l, p = 0, q =0, len = 0, left = 0, right = 0;
    
    for (i  = 0; i < alen - len; i++)
    {
        for (j = 0; j < blen - len; j++)
        {
            if (a[i] == b[j] && a[i + len] == b[j + len])
            {
                for (k = i + 1, l = j + 1; 
                     a[k] == b[l] && k < alen && l < blen
                     ; k++, l++);

                if (k - i > len)
                {
                    p = i;
                    q = j;
                    len = k - i;
                } 
            }
        }
    }
    
    if (!len) return 0;

    if (p !=0 && q != 0)
        left = rsimil(a, p, b, q);
    
    i = p + len;
    alen -= i;
    j = q + len;
    blen  -= j;
    
    if (alen != 0 && blen != 0)
        right = rsimil(a + i, alen, b + j, blen);

    return len + left + right;
}

int simil(const char* a, const char* b)
{
    int alen = std::strlen(a);
    int blen = std::strlen(b);

    if (!alen || !blen) return 0;

    return (rsimil(a, alen, b, blen) * 200) / (alen + blen);
}

int similarity_score(const std::string& a, const std::string& b)
{
    return simil(a.c_str(), b.c_str());
}
