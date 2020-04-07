

#ifdef WIN32
#include <hash_map>
using std::hash_map;
#else
#include <ext/hash_map>
using __gnu_cxx::hash_map;
#endif
