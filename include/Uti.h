
// copy the string and return a char pointer
inline void split(const string &s, const char *delim, vector <string> &v) {
  char *dup = strdup(s.c_str());
  char *token = strtok(dup, delim);
  while (token != NULL) {
    v.push_back(string(token));
    // the call is treated as a subsequent calls to strtok:
    // the function continues from where it left in previous invocation
    token = strtok(NULL, delim);
  }
  free(dup);
}