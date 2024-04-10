#include <utils.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

void *memdup(const void *src, const size_t n) {
  void *dest = malloc(n);
  if(dest == NULL)
    return NULL;
  return memcpy(dest,src,n);
}

int int_compar(const void *a, const void *b) {
  return int_compar_r(a, b, NULL);
}

int int_compar_r(const void *a, const void *b, void *arg) {
  int order = (arg == NULL) ? 1 : *((const int *) arg);
  if(*((int *) a) < *((int *) b))
    return -order;
  if(*((int *) a) > *((int *) b))
    return order;
  return 0;
}

int dbl_compar(const void *a, const void *b) {
  return dbl_compar_r(a, b, NULL);
}

int dbl_compar_r(const void *a, const void *b, void *arg) {
  int order = (arg == NULL) ? 1 : *((const int *) arg);
  if(*((double *) a) < *((double *) b))
    return -order;
  if(*((double *) a) > *((double *) b))
    return order;
  return 0;
}

const void *min(const void *_x, const size_t nmemb, const size_t size, int (*compar)(const void *, const void *)) {
  const char (*x)[size] = (const char (*)[size]) _x;
  const void *res = _x;

  if(nmemb == 0) return NULL;
  for(size_t i = 1; i < nmemb; i++) {
    if(compar(&x[i], res) < 0)
      res = &x[i];
  }

  return res;
}

const void *max(const void *_x, const size_t nmemb, const size_t size, int (*compar)(const void *, const void *)) {
  const char (*x)[size] = (const char (*)[size]) _x;
  const void *res = _x;

  if(nmemb == 0) return NULL;
  for(size_t i = 1; i < nmemb; i++) {
    if(compar(&x[i], res) > 0)
      res = &x[i];
  }

  return res;
}

const void *min_r(const void *_x, const size_t nmemb, const size_t size, int (*compar)(const void *, const void *, void *), void *arg) {
  const char (*x)[size] = (const char (*)[size]) _x;
  const void *res = _x;

  if(nmemb == 0) return NULL;
  for(size_t i = 1; i < nmemb; i++) {
    if(compar(&x[i], res, arg) < 0)
      res = &x[i];
  }

  return res;
}

const void *max_r(const void *_x, const size_t nmemb, const size_t size, int (*compar)(const void *, const void *, void *), void *arg) {
  const char (*x)[size] = (const char (*)[size]) _x;
  const void *res = _x;

  if(nmemb == 0) return NULL;
  for(size_t i = 1; i < nmemb; i++) {
    if(compar(&x[i], res, arg) > 0)
      res = &x[i];
  }

  return res;
}

struct _argsort_arg_t {
  const void *base;
  const size_t base_size;
  int (*base_compar)(const void *, const void *);
};

int _argsort_compar(const void *a, const void *b, void *arg) {
  const struct _argsort_arg_t *c = (const struct _argsort_arg_t *) arg;
  return c->base_compar(c->base+*((size_t *)a)*c->base_size,
                        c->base+*((size_t *)b)*c->base_size);
}

void argsort(size_t *idx, const void *base, const size_t nmemb, const size_t size, int (*compar)(const void *, const void *)) {
  struct _argsort_arg_t c = {base, size, compar};
  for(size_t i = 0; i < nmemb; i++) idx[i] = i;
  qsort_r(idx, nmemb, sizeof(size_t), _argsort_compar, &c);
}

struct _argsort_r_arg_t {
  const void *base;
  const size_t base_size;
  int (*base_compar)(const void *, const void *, void *);
  void *base_arg;
};

int _argsort_r_compar(const void *a, const void *b, void *arg) {
  const struct _argsort_r_arg_t *c = (const struct _argsort_r_arg_t *) arg;
  return c->base_compar(c->base+*((size_t *)a)*c->base_size,
                        c->base+*((size_t *)b)*c->base_size,
                        c->base_arg);
}

void argsort_r(size_t *idx, const void *base, const size_t nmemb, const size_t size, int (*compar)(const void *, const void *, void *), void *arg) {
  struct _argsort_r_arg_t c = {base, size, compar, arg};
  for(size_t i = 0; i < nmemb; i++) idx[i] = i;
  qsort_r(idx, nmemb, sizeof(size_t), _argsort_r_compar, &c);
}

int mrand_r (unsigned int *seed) {
  // Glibc 2.24 implementation of rand_r
  unsigned int next = *seed;
  int result;

  next *= 1103515245;
  next += 12345;
  result = (unsigned int) (next / 65536) % 2048;

  next *= 1103515245;
  next += 12345;
  result <<= 10;
  result ^= (unsigned int) (next / 65536) % 1024;

  next *= 1103515245;
  next += 12345;
  result <<= 10;
  result ^= (unsigned int) (next / 65536) % 1024;

  *seed = next;

  return result;
}

double mrandn(unsigned int *seed) {
  double U1, U2, W, s, X1;
  static double X2;
  static char call = 0;

  if (call) {
    call = 0;
    return X2;
  }

  do {
    U1 = 2. * mrand_r(seed) / MRAND_MAX - 1;
    U2 = 2. * mrand_r(seed) / MRAND_MAX - 1;
    W = U1 * U1 + U2 * U2;
  } while (W == 0 || 1 <= W);

  s = sqrt(-2. * log (W) / W);
  X1 = U1 * s;
  X2 = U2 * s;

  call = 1;

  return X1;
}

static const unsigned int _crc32_table[256] = {
  0x00000000U,0x04C11DB7U,0x09823B6EU,0x0D4326D9U,0x130476DCU,0x17C56B6BU,0x1A864DB2U,0x1E475005U,
  0x2608EDB8U,0x22C9F00FU,0x2F8AD6D6U,0x2B4BCB61U,0x350C9B64U,0x31CD86D3U,0x3C8EA00AU,0x384FBDBDU,
  0x4C11DB70U,0x48D0C6C7U,0x4593E01EU,0x4152FDA9U,0x5F15ADACU,0x5BD4B01BU,0x569796C2U,0x52568B75U,
  0x6A1936C8U,0x6ED82B7FU,0x639B0DA6U,0x675A1011U,0x791D4014U,0x7DDC5DA3U,0x709F7B7AU,0x745E66CDU,
  0x9823B6E0U,0x9CE2AB57U,0x91A18D8EU,0x95609039U,0x8B27C03CU,0x8FE6DD8BU,0x82A5FB52U,0x8664E6E5U,
  0xBE2B5B58U,0xBAEA46EFU,0xB7A96036U,0xB3687D81U,0xAD2F2D84U,0xA9EE3033U,0xA4AD16EAU,0xA06C0B5DU,
  0xD4326D90U,0xD0F37027U,0xDDB056FEU,0xD9714B49U,0xC7361B4CU,0xC3F706FBU,0xCEB42022U,0xCA753D95U,
  0xF23A8028U,0xF6FB9D9FU,0xFBB8BB46U,0xFF79A6F1U,0xE13EF6F4U,0xE5FFEB43U,0xE8BCCD9AU,0xEC7DD02DU,
  0x34867077U,0x30476DC0U,0x3D044B19U,0x39C556AEU,0x278206ABU,0x23431B1CU,0x2E003DC5U,0x2AC12072U,
  0x128E9DCFU,0x164F8078U,0x1B0CA6A1U,0x1FCDBB16U,0x018AEB13U,0x054BF6A4U,0x0808D07DU,0x0CC9CDCAU,
  0x7897AB07U,0x7C56B6B0U,0x71159069U,0x75D48DDEU,0x6B93DDDBU,0x6F52C06CU,0x6211E6B5U,0x66D0FB02U,
  0x5E9F46BFU,0x5A5E5B08U,0x571D7DD1U,0x53DC6066U,0x4D9B3063U,0x495A2DD4U,0x44190B0DU,0x40D816BAU,
  0xACA5C697U,0xA864DB20U,0xA527FDF9U,0xA1E6E04EU,0xBFA1B04BU,0xBB60ADFCU,0xB6238B25U,0xB2E29692U,
  0x8AAD2B2FU,0x8E6C3698U,0x832F1041U,0x87EE0DF6U,0x99A95DF3U,0x9D684044U,0x902B669DU,0x94EA7B2AU,
  0xE0B41DE7U,0xE4750050U,0xE9362689U,0xEDF73B3EU,0xF3B06B3BU,0xF771768CU,0xFA325055U,0xFEF34DE2U,
  0xC6BCF05FU,0xC27DEDE8U,0xCF3ECB31U,0xCBFFD686U,0xD5B88683U,0xD1799B34U,0xDC3ABDEDU,0xD8FBA05AU,
  0x690CE0EEU,0x6DCDFD59U,0x608EDB80U,0x644FC637U,0x7A089632U,0x7EC98B85U,0x738AAD5CU,0x774BB0EBU,
  0x4F040D56U,0x4BC510E1U,0x46863638U,0x42472B8FU,0x5C007B8AU,0x58C1663DU,0x558240E4U,0x51435D53U,
  0x251D3B9EU,0x21DC2629U,0x2C9F00F0U,0x285E1D47U,0x36194D42U,0x32D850F5U,0x3F9B762CU,0x3B5A6B9BU,
  0x0315D626U,0x07D4CB91U,0x0A97ED48U,0x0E56F0FFU,0x1011A0FAU,0x14D0BD4DU,0x19939B94U,0x1D528623U,
  0xF12F560EU,0xF5EE4BB9U,0xF8AD6D60U,0xFC6C70D7U,0xE22B20D2U,0xE6EA3D65U,0xEBA91BBCU,0xEF68060BU,
  0xD727BBB6U,0xD3E6A601U,0xDEA580D8U,0xDA649D6FU,0xC423CD6AU,0xC0E2D0DDU,0xCDA1F604U,0xC960EBB3U,
  0xBD3E8D7EU,0xB9FF90C9U,0xB4BCB610U,0xB07DABA7U,0xAE3AFBA2U,0xAAFBE615U,0xA7B8C0CCU,0xA379DD7BU,
  0x9B3660C6U,0x9FF77D71U,0x92B45BA8U,0x9675461FU,0x8832161AU,0x8CF30BADU,0x81B02D74U,0x857130C3U,
  0x5D8A9099U,0x594B8D2EU,0x5408ABF7U,0x50C9B640U,0x4E8EE645U,0x4A4FFBF2U,0x470CDD2BU,0x43CDC09CU,
  0x7B827D21U,0x7F436096U,0x7200464FU,0x76C15BF8U,0x68860BFDU,0x6C47164AU,0x61043093U,0x65C52D24U,
  0x119B4BE9U,0x155A565EU,0x18197087U,0x1CD86D30U,0x029F3D35U,0x065E2082U,0x0B1D065BU,0x0FDC1BECU,
  0x3793A651U,0x3352BBE6U,0x3E119D3FU,0x3AD08088U,0x2497D08DU,0x2056CD3AU,0x2D15EBE3U,0x29D4F654U,
  0xC5A92679U,0xC1683BCEU,0xCC2B1D17U,0xC8EA00A0U,0xD6AD50A5U,0xD26C4D12U,0xDF2F6BCBU,0xDBEE767CU,
  0xE3A1CBC1U,0xE760D676U,0xEA23F0AFU,0xEEE2ED18U,0xF0A5BD1DU,0xF464A0AAU,0xF9278673U,0xFDE69BC4U,
  0x89B8FD09U,0x8D79E0BEU,0x803AC667U,0x84FBDBD0U,0x9ABC8BD5U,0x9E7D9662U,0x933EB0BBU,0x97FFAD0CU,
  0xAFB010B1U,0xAB710D06U,0xA6322BDFU,0xA2F33668U,0xBCB4666DU,0xB8757BDAU,0xB5365D03U,0xB1F740B4U,
};

unsigned int hashcode(const void *data, const size_t n) {
  const unsigned char *cdata = (const unsigned char *) data;
  unsigned int h = 0;

  for(size_t i = 0; i < n; i++)
    h = _crc32_table[cdata[i] ^ ((h >> 24) & 0xFF)] ^ (h << 8);

  return h;
}

void swap(void *a, void *b, size_t size) {
  uint8_t *_a = a, *_b = b, tmp;
  size_t i;

  for(i = 0; i < size; i++) {
    tmp = _a[i];
    _a[i] = _b[i];
    _b[i] = tmp;
  }
}

void shuffle(void *base, size_t nmemb, size_t size, unsigned int *seed) {
  uint8_t (*a)[size] = (uint8_t (*)[size]) base;
  size_t i, j;

  for(i = nmemb - 1; 0 < i; i--){
    j = mrand_r(seed) % (i + 1);
    swap(a[i], a[j], size);
  }
}

int parse_double(double *d, const char *s, char **endptr, int last) {
  char *_endptr;

  *d = strtod(s, &_endptr);
  if(s == _endptr || (last && *_endptr != '\0')) return -1;
  if(endptr) *endptr = _endptr;

  return 0;
}
