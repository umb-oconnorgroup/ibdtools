static inline int
strcmp_i(const char *a, const char *b)
{
    while (1) {
        if (*a != *b)
            return *a - *b;
        else {
            if (*a == '\0')
                return 0;
            ++a;
            ++b;
        }
    }
}
