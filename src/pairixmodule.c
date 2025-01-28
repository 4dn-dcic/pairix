/*-
 * The MIT License
 *
 * Copyright (c) 2011 Seoul National University.
 *               2016 Soo Lee
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/*
 * Contact: Soo Lee (duplexa@gmail.com)
 */

#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include "pairix.h"

static PyObject *PairixError;
static PyObject *PairixWarning;

typedef struct {
    PyObject_HEAD
    pairix_t *tb;
    char *fn;
    int linecount;
} PairixObject;

typedef struct {
    PyObject_HEAD
    PairixObject *tbobj;
    sequential_iter_t *iter;
} PairixIteratorObject;


static PyTypeObject Pairix_Type, PairixIterator_Type;

/* --- PairixIterator --------------------------------------------------- */

static PyObject *
pairixiter_create(PairixObject *parentidx, sequential_iter_t *iter)
{
    PairixIteratorObject *self;

    self = PyObject_New(PairixIteratorObject, &PairixIterator_Type);
    if (self == NULL)
        return NULL;

    Py_INCREF(parentidx);
    self->tbobj = parentidx;
    self->iter = iter;

    return (PyObject *)self;
}

static void
pairixiter_dealloc(PairixIteratorObject *self)
{
    destroy_sequential_iter(self->iter);
    Py_DECREF(self->tbobj);
    PyObject_Del(self);
}

static PyObject *
pairixiter_iter(PyObject *self)
{
    Py_INCREF(self);
    return self;
}

#if PY_MAJOR_VERSION < 3
# define PYOBJECT_FROM_STRING_AND_SIZE PyString_FromStringAndSize
#else
# define PYOBJECT_FROM_STRING_AND_SIZE PyUnicode_FromStringAndSize
#endif

#if PY_MAJOR_VERSION < 3
# define PYOBJECT_FROM_STRING PyString_FromString
#else
# define PYOBJECT_FROM_STRING PyUnicode_FromString
#endif

static PyObject *
pairixiter_iternext(PairixIteratorObject *self)
{
    const char *chunk;
    int len, i;
    char delimiter = ti_get_delimiter(self->tbobj->tb->idx);

    chunk = sequential_ti_read(self->iter, &len);

    if (chunk != NULL) {
        PyObject *ret, *column;
        Py_ssize_t colidx;
        const char *ptr, *begin;

        ret = PyList_New(0);
        if (ret == NULL)
            return NULL;

        colidx = 0;
        ptr = begin = chunk;
        for (i = len; i > 0; i--, ptr++)
            if (*ptr == delimiter) {
                column = PYOBJECT_FROM_STRING_AND_SIZE(begin,
                                                       (Py_ssize_t)(ptr - begin));
                if (column == NULL || PyList_Append(ret, column) == -1) {
                    Py_DECREF(ret);
                    return NULL;
                }

                Py_DECREF(column);
                begin = ptr + 1;
                colidx++;
            }

        column = PYOBJECT_FROM_STRING_AND_SIZE(begin, (Py_ssize_t)(ptr - begin));
        if (column == NULL || PyList_Append(ret, column) == -1) {
            Py_DECREF(ret);
            return NULL;
        }
        Py_DECREF(column);

        return ret;
    }
    else
        return NULL;
}

static PyMethodDef pairixiter_methods[] = {
    {NULL, NULL} /* sentinel */
};

static PyTypeObject PairixIterator_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "pypairix.iter",      /*tp_name*/
    sizeof(PairixIteratorObject), /*tp_basicsize*/
    0,                          /*tp_itemsize*/
    /* methods */
    (destructor)pairixiter_dealloc,  /*tp_dealloc*/
    0,                          /*tp_print*/
    0,                          /*tp_getattr*/
    0,                          /*tp_setattr*/
    0,                          /*tp_compare*/
    0,                          /*tp_repr*/
    0,                          /*tp_as_number*/
    0,                          /*tp_as_sequence*/
    0,                          /*tp_as_mapping*/
    0,                          /*tp_hash*/
    0,                          /*tp_call*/
    0,                          /*tp_str*/
    0,                          /*tp_getattro*/
    0,                          /*tp_setattro*/
    0,                          /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,         /*tp_flags*/
    0,                          /*tp_doc*/
    0,                          /*tp_traverse*/
    0,                          /*tp_clear*/
    0,                          /*tp_richcompare*/
    0,                          /*tp_weaklistoffset*/
    pairixiter_iter,             /*tp_iter*/
    (iternextfunc)pairixiter_iternext, /*tp_iternext*/
    pairixiter_methods,          /*tp_methods*/
    0,                          /*tp_members*/
    0,                          /*tp_getset*/
    0,                          /*tp_base*/
    0,                          /*tp_dict*/
    0,                          /*tp_descr_get*/
    0,                          /*tp_descr_set*/
    0,                          /*tp_dictoffset*/
    0,                          /*tp_init*/
    0,                          /*tp_alloc*/
    0,                          /*tp_new*/
    0,                          /*tp_free*/
    0,                          /*tp_is_gc*/
};

/* ------------------------build_index------------------------- */

int build_index(char *inputfilename, char *preset, int sc, int bc, int ec, int sc2, int bc2, int ec2, char *delimiter, char *meta_char, int line_skip, char *region_split_character, int force, int zero){

  if(force==0){
    char *fnidx = calloc(strlen(inputfilename) + 5, 1);
    strcat(strcpy(fnidx, inputfilename), ".px2");
    struct stat stat_px2, stat_input;
     if (stat(fnidx, &stat_px2) == 0){  // file exists.
         // Before complaining, check if the input file isn't newer. This is a common source of errors,
         // people tend not to notice that pairix failed
         stat(inputfilename, &stat_input);
         if ( stat_input.st_mtime <= stat_px2.st_mtime ) { free(fnidx); return(-4);}
     }
     free(fnidx);
  }
  if ( bgzf_is_bgzf(inputfilename)!=1 ) return(-3);
  else {
    ti_conf_t conf;
    if (strcmp(preset, "") == 0 && sc == 0 && bc == 0){
      int l = strlen(inputfilename);
      int strcasecmp(const char *s1, const char *s2);
      if (l>=7 && strcasecmp(inputfilename+l-7, ".gff.gz") == 0) conf = ti_conf_gff;
      else if (l>=7 && strcasecmp(inputfilename+l-7, ".bed.gz") == 0) conf = ti_conf_bed;
      else if (l>=7 && strcasecmp(inputfilename+l-7, ".sam.gz") == 0) conf = ti_conf_sam;
      else if (l>=7 && strcasecmp(inputfilename+l-7, ".vcf.gz") == 0) conf = ti_conf_vcf;
      else if (l>=10 && strcasecmp(inputfilename+l-10, ".psltbl.gz") == 0) conf = ti_conf_psltbl;
      else if (l>=9 && strcasecmp(inputfilename+l-9, ".pairs.gz") == 0) conf = ti_conf_pairs;
      else return(-5); // file extension not recognized and no preset specified
    }
    else if (strcmp(preset, "") == 0 && sc != 0 && bc != 0){
      conf.sc = sc;
      conf.bc = bc;
      conf.ec = ec;
      conf.sc2 = sc2;
      conf.bc2 = bc2;
      conf.ec2 = ec2;
      conf.delimiter = delimiter[0];
      conf.region_split_character = region_split_character[0];
      conf.meta_char = (int)(meta_char[0]);
      conf.line_skip = line_skip;
    }
    else if (strcmp(preset, "gff") == 0) conf = ti_conf_gff;
    else if (strcmp(preset, "bed") == 0) conf = ti_conf_bed;
    else if (strcmp(preset, "sam") == 0) conf = ti_conf_sam;
    else if (strcmp(preset, "vcf") == 0 || strcmp(preset, "vcf4") == 0) conf = ti_conf_vcf;
    else if (strcmp(preset, "psltbl") == 0) conf = ti_conf_psltbl;
    else if (strcmp(preset, "pairs") == 0) conf = ti_conf_pairs;
    else if (strcmp(preset, "merged_nodups") == 0) conf = ti_conf_merged_nodups;
    else if (strcmp(preset, "old_merged_nodups") == 0) conf = ti_conf_old_merged_nodups;
    else return(-2);  // wrong preset

    if(sc!=0) conf.sc=sc;
    if(bc!=0) conf.bc=bc;
    if(ec!=0) conf.ec=ec;
    if(sc2!=0) conf.sc2=sc2;
    if(bc2!=0) conf.bc2=bc2;
    if(ec2!=0) conf.ec2=ec2;
    if(line_skip!=0) conf.line_skip=line_skip;
    if(strcmp(delimiter, "\t")) conf.delimiter = delimiter[0];
    if(strcmp(meta_char, "#")) conf.meta_char = meta_char[0];
    if(strcmp(region_split_character, DEFAULT_REGION_SPLIT_CHARACTER_STR)) conf.region_split_character = region_split_character[0];


    if(zero) conf.preset |= TI_FLAG_UCSC; // zero-based indexing

    return ti_index_build(inputfilename, &conf);  // -1 if failed
  }
}

static char indexer_docstring[] = (
"build_index(filename, preset, sc, bc, ec, sc2, bc2, ec2, delimiter, meta_char, region_split_character, force)\n\n"
"Build pairix index for a pairs file.\n\n"
"Parameters\n"
"----------\n"
"filename : str\n"
"    Name of bgzipped file to index.\n"
"preset : str, optional (default : '')\n"
"    One of the following strings: 'gff', 'bed', 'sam', 'vcf', 'psltbl' \n"
"   (1D-indexing) or 'pairs', 'merged_nodups', 'old_merged_nodups' (2D-indexing).\n"
"sc, bc, ec : int, optional\n"
"    First sequence, start, and end column indexes (1-based). Zero (0) means \n"
"    not specified. If preset is given, preset overrides these. If preset is \n"
"    not given, sc2 and bc2 are required.\n"
"sc2, bc2, ec2 : int, optional\n"
"    Second sequence, start, and end column indexes (1-based). Zero (0) \n"
"    means not specified. If preset is given, preset overrides these. If sc, \n"
"    bc are specified but not sc2 and bc2, the file will be 1D-indexed.\n"
"delimiter : str, optional (default: '\\t')\n"
"    Data delimiter, e.g. '\\t' or ' '. If preset is given, preset overrides \n"
"    delimiter.\n"
"meta_char : str, optional (default: '#')\n"
"    Comment character. Lines beginning with this character are skipped when \n"
"    creating an index. If preset is given, preset overrides comment_char.\n"
"line_skip : int, optional (default: 0)\n"
"    Number of lines to skip in the beginning\n"
"region_split_character : char, option (default: '|')\n"
"    Character used to separate two regions.\n"
"force : int, optional (default: 0)\n"
"    If 1, overwrite existing index file. If 0, do not overwrite unless the \n"
"    index file is older than the bgzipped file.\n"
"zero : int, optional (default 0)\n"
"    If 1, create a zero-based index. (default one-based)\n");


static PyObject *indexer_build_index(PyObject *self, PyObject *args, PyObject *kwargs)
{
    // args: char *inputfilename, char *preset, int sc, int bc, int ec, int sc2, int bc2, int ec2, char *delimiter, char *meta_char, int line_skip, char *region_split_character, int force, int zero
    char *inputfilename, *preset, *delimiter, *meta_char, *region_split_character;
    int sc, bc, ec, sc2, bc2, ec2, line_skip, force, zero, result;
    // default arg values
    preset="";
    delimiter="\t";
    meta_char="#";
    region_split_character=DEFAULT_REGION_SPLIT_CHARACTER_STR;
    sc = 0;
    bc = 0;
    ec = 0;
    sc2 = 0;
    bc2 = 0;
    ec2 = 0;
    line_skip = 0;
    force = 0;
    zero = 0;

    static char *kwlist[] = {"inputfilename", "preset", "sc", "bc", "ec", "sc2", "bc2", "ec2", "delimiter", "meta_char", "line_skip", "region_split_character", "force", "zero", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s|siiiiiissisii", kwlist, &inputfilename, &preset, &sc, &bc, &ec, &sc2, &bc2, &ec2, &delimiter, &meta_char, &line_skip, &region_split_character, &force, &zero)){
        PyErr_SetString(PairixError, "Argument error! build_index() requires the following args:\n<filename (str)>.\nOptional args:\n<preset (str)> one of the following strings: 'gff', 'bed', 'sam', 'vcf', 'psltbl' (1D-indexing) or 'pairs', 'merged_nodups', 'old_merged_nodups' (2D-indexing). If preset is '', at least some of the custom parameters must be given instead (sc, bc, ec, sc2, bc2, ec2, delimiter, comment_char, line_skip, region_split_character). (default '')\n<sc (int)> first sequence (chromosome) column index (1-based). Zero (0) means not specified. If preset is given, preset overrides sc. If preset is not given, this one is required. (default 0)\n<bc (int)> first start position column index (1-based). Zero (0) means not specified. If preset is given, preset overrides bc. If preset is not given, this one is required. (default 0)\n<ec (int)> first end position column index (1-based). Zero (0) means not specified. If preset is given, preset overrides ec. (default 0)\n<sc2 (int)> second sequence (chromosome) column index (1-based). Zero (0) means not specified. If preset is given, preset overrides sc2. If sc, bc are specified but not sc2 and bc2, it is 1D-indexed. (default 0)\n<bc2 (int)> second start position column index (1-based). Zero (0) means not specified. If preset is given, preset overrides bc2. (default 0)\n<ec2 (int)> second end position column index (1-based). Zero (0) means not specified. If preset is given, preset overrides ec2. (default 0)\n<delimiter (str)> delimiter (e.g. '\\t' or ' ') (default '\\t'). If preset is given, preset overrides delimiter.\n<meta_char (str)> comment character. Lines beginning with this character are skipped when creating an index. If preset is given, preset overrides comment_char (default '#')\n<line_skip (int)> number of lines to skip in the beginning. (default 0)\n<region_split_character (char)> character used to separate two regions. (default '|')\n<force (int)> If 1, overwrite existing index file. If 0, do not overwrite unless the index file is older than the bgzipped file. (default 0). <zero (int)> If 1, create a zero-based index. (default one-based)\n");
        return NULL;
    }

    result = build_index(inputfilename, preset, sc, bc, ec, sc2, bc2, ec2, delimiter, meta_char, line_skip, region_split_character, force, zero);

    if (result == -1) {
        PyErr_SetString(PairixError, "Can't create index.");
        return NULL;
    }
    else if (result == -2) {
        PyErr_SetString(PairixError, "Can't recognize preset.");
        return NULL;
    }
    else if (result == -3) {
        PyErr_SetString(PairixError, "Was bgzip used to compress this file?");
        return NULL;
    }
    else if (result == -4) {
        PyErr_SetString(PairixError, "The index file exists. Please use force=1 to overwrite.");
        return NULL;
    }
    else if (result == -5) {
        PyErr_SetString(PairixError, "Can't recognize file type, with no preset specified.");
        return NULL;
    }

    Py_RETURN_NONE;
}

static PyMethodDef indexer_methods[] = {
    {"build_index", (PyCFunction)indexer_build_index, METH_VARARGS | METH_KEYWORDS, indexer_docstring},
    {NULL, NULL} /* sentinel */
};

/* --- Pairix ----------------------------------------------------------- */

static PyObject *
pairix_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PairixObject *self;
    const char *fn, *fnidx=NULL;
    static char *kwnames[]={"fn", "fnidx", NULL};
    pairix_t *tb;
    int i;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|z:open",
                                     kwnames, &fn, &fnidx))
        return NULL;

    tb = ti_open(fn, fnidx);
    if (tb == NULL) {
        PyErr_SetString(PairixError, "Can't open the index file.");
        return NULL;
    }

    self = (PairixObject *)type->tp_alloc(type, 0);
    if (self == NULL)
        return NULL;

    self->tb = tb;
    self->fn = strdup(fn);
    self->tb->idx = ti_index_load(self->fn);
    if(self->tb->idx == NULL) {
        PyErr_SetString(PairixError, "Can't open the index file.");
        return NULL;
    }
    self->linecount = get_linecount(self->tb->idx);
    return (PyObject *)self;
}

static void
pairix_dealloc(PairixObject *self)
{
    free(self->fn);
    ti_close(self->tb);
    PyObject_Del(self);
}


static PyObject *
pairix_query(PairixObject *self, PyObject *args)
{
    char *name;
    int begin, end, tid_test;
    sequential_iter_t *result;

    if (!PyArg_ParseTuple(args, "sii:query", &name, &begin, &end)){
        PyErr_SetString(PairixError, "Argument error! query() takes the following args: <chromosome (str)> <begin (int)> <end (int)>");
        return pairixiter_create(self, NULL);
    }

    tid_test = ti_query_tid(self->tb, name, begin, end);
    if (tid_test == -1) {
        PyErr_WarnEx(PairixWarning, "Cannot find matching chromosome name. Check that chromosome naming conventions match between your query and input file.",1);
        return pairixiter_create(self, NULL);
    }
    else if (tid_test == -2){
        PyErr_WarnEx(PairixWarning, "The start coordinate must be less than the end coordinate.", 1);
        return pairixiter_create(self, NULL);
    }
    else if (tid_test == -3){
        PyErr_WarnEx(PairixWarning, "The specific cause could not be found. Please adjust your arguments.", 1);
        return pairixiter_create(self, NULL);
    }

    result = ti_query_general(self->tb, name, begin, end);
    return pairixiter_create(self, result); // result may be null but that's okay.
}

static PyObject *
pairix_queryi(PairixObject *self, PyObject *args)
{
    int tid, begin, end;
    sequential_iter_t *result;

    if (!PyArg_ParseTuple(args, "iii:queryi", &tid, &begin, &end)){
        PyErr_SetString(PairixError, "Argument error! queryi() takes three integers: tid, begin and end");
        return pairixiter_create(self, NULL);
    }

    result = ti_queryi_general(self->tb, tid, begin, end);
    return pairixiter_create(self, result);  // result may be null but that's okay
}

static PyObject *
pairix_querys(PairixObject *self, PyObject *args)
{
    const char *reg;
    int tid_test;
    sequential_iter_t *result;

    if (!PyArg_ParseTuple(args, "s:querys", &reg)){
        PyErr_SetString(PairixError, "Argument error! querys2D() takes one str formatted as: '{CHR}:{START}-{END}'");
        return pairixiter_create(self, NULL);
    }

    tid_test = ti_querys_tid(self->tb, reg);
    if (tid_test == -1) {
        PyErr_WarnEx(PairixWarning, "Cannot find matching chromosome name. Check that chromosome naming conventions match between your query and input file.", 1);
        return pairixiter_create(self, NULL);
    }
    else if (tid_test == -2){
        PyErr_WarnEx(PairixWarning, "The start coordinate must be less than the end coordinate.", 1);
        return pairixiter_create(self, NULL);
    }
    else if (tid_test == -3){
        PyErr_WarnEx(PairixWarning, "The specific cause could not be found. Please adjust your arguments.", 1);
        return pairixiter_create(self, NULL);
    }

    result = ti_querys_general(self->tb, reg);
    return pairixiter_create(self, result); // result may be null but that's okay
}

/* ------- PAIRIX 2D QUERYING METHODS ------- */

static PyObject *
pairix_query_2D(PairixObject *self, PyObject *args)
{
    char *name, *name2;
    int begin, end, begin2, end2, tid_test, tid_test_rev, flip;
    sequential_iter_t *result;
    flip = 0;

    if (!PyArg_ParseTuple(args, "siisii|i:query2D", &name, &begin, &end, &name2, &begin2, &end2, &flip)){
        PyErr_SetString(PairixError, "Argument error! query2D() takes the following args: <1st_chromosome (str)> <begin (int)> <end (int)> <2nd_chromosome (str)> <begin (int)> <end (int)> [<autoflip>]. Optionally, include an integer = 1 to test chromsomes in flipped order on an error (autoflip).");
        return pairixiter_create(self, NULL);
    }

    tid_test = ti_query_2d_tid(self->tb, name, begin, end, name2, begin2, end2);
    if (tid_test == -1) {
        // Test if reversing chromosome order fixes the issues
        tid_test_rev = ti_query_2d_tid(self->tb, name2, begin2, end2, name, begin, end);
        if (tid_test_rev != -1 && tid_test_rev != -2 && tid_test_rev != -3) {
            result = ti_query_2d_general(self->tb, name2, begin2, end2, name, begin, end);
            if (flip == 1){
                if (result == NULL) {
                   // PyErr_SetString(PairixError, "Input error! Cannot find matching chromosome names. Check that chromosome naming conventions match between your query and input file.");
                   PyErr_WarnEx(PairixWarning, "Cannot find matching chromosome pair. Check that chromosome naming conventions match between your query and input file.",1);
                   return pairixiter_create(self, NULL);
                }else{
                    return pairixiter_create(self, result);
                }
            }
            else{
                // PyErr_SetString(PairixError, "Input error! Cannot find matching chromosome names. Check that chromosome naming conventions match between your query and input file. You may wish to also automatically test chromsomes in flipped order. To do this, include 1 as the last argument.");
                PyErr_WarnEx(PairixWarning, "Cannot find matching chromosome pair. Check that chromosome naming conventions match between your query and input file. You may wish to also automatically test chromsomes in flipped order. To do this, include 1 as the last argument.",1);
                return pairixiter_create(self, NULL);
            }
        }
    }
    else if (tid_test == -2){
        PyErr_WarnEx(PairixWarning, "The start coordinate must be less than the end coordinate.", 1);
        return pairixiter_create(self, NULL);
    }
    else if (tid_test == -3){
        PyErr_WarnEx(PairixWarning, "The specific cause could not be found. Please adjust your arguments.", 1);
        return pairixiter_create(self, NULL);
    }

    result = ti_query_2d_general(self->tb, name, begin, end, name2, begin2, end2);
    return pairixiter_create(self, result); // result may be null but that's okay
}

static PyObject *
pairix_queryi_2D(PairixObject *self, PyObject *args)
{
    int tid, begin, end, begin2, end2;
    sequential_iter_t *result;

    if (!PyArg_ParseTuple(args, "iiiii:queryi2D", &tid, &begin, &end, &begin2, &end2)){
        PyErr_SetString(PairixError, "Argument error! queryi2D() takes five integers: tid, begin, end, begin2 and end2");
        return pairixiter_create(self, NULL);
    }

    result = ti_queryi_2d_general(self->tb, tid, begin, end, begin2, end2);
    return pairixiter_create(self, result);  // result may be null but that's okay
}

static PyObject *
pairix_querys_2D(PairixObject *self, PyObject *args)
{
    const char *reg, *reg2;
    int tid_test, tid_test_rev, flip;
    sequential_iter_t *result;

    flip = 0;

    if (!PyArg_ParseTuple(args, "s|i:querys2D", &reg, &flip)){
        PyErr_SetString(PairixError, "Argument error! querys2D() takes the following args: <query_str> [<autoflip>]. Query_str is a str formatted as: '{CHR}:{START}-{END}|{CHR}:{START}-{END}' Optionally, include an integer = 1 to test chromsomes in flipped order on an error (autoflip).");
        return pairixiter_create(self, NULL);
    }

    tid_test = ti_querys_tid(self->tb, reg);
    if (tid_test == -1) {
        reg2 = flip_region(reg, get_region_split_character(self->tb));
        tid_test_rev = ti_querys_tid(self->tb, reg2);
        if (tid_test_rev != -1 && tid_test_rev != -2 && tid_test_rev != -3) {
            result = ti_querys_2d_general(self->tb, reg2);
            if (flip == 1){
                if (result == NULL) {
                   // PyErr_SetString(PairixError, "Input error! Cannot find matching chromosome names. Check that chromosome naming conventions match between your query and input file.");
                   PyErr_WarnEx(PairixWarning, "Cannot find matching chromosome pair. Check that chromosome naming conventions match between your query and input file.",1);
                   return pairixiter_create(self, NULL);
                }else{
                    return pairixiter_create(self, result);
                }
            }
            else{
                // PyErr_SetString(PairixError, "Input error! Cannot find matching chromosome names. Check that chromosome naming conventions match between your query and input file. You may wish to also automatically test chromsomes in flipped order. To do this, include 1 as the last argument.");
                PyErr_WarnEx(PairixWarning, "Cannot find matching chromosome pair. Check that chromosome naming conventions match between your query and input file. You may wish to also automatically test chromsomes in flipped order. To do this, include 1 as the last argument.",1);
                return pairixiter_create(self, NULL);
            }
        }
    }
    else if (tid_test == -2){
        PyErr_WarnEx(PairixWarning, "The start coordinate must be less than the end coordinate.", 1);
        return pairixiter_create(self, NULL);
    }
    else if (tid_test == -3){
        PyErr_WarnEx(PairixWarning, "The specific cause could not be found. Please adjust your arguments.", 1);
        return pairixiter_create(self, NULL);
    }

    result = ti_querys_2d_general(self->tb, reg);
    return pairixiter_create(self, result);  // result may be null but that's okay
}


static PyObject *
pairix_get_linecount(PairixObject *self)
{
  return Py_BuildValue("l", self->linecount);
}

static PyObject *
pairix_get_blocknames(PairixObject *self)
{
  int n,i;
  char **blocknames = (char**)ti_seqname(self->tb->idx, &n);
  PyObject *bnames = PyList_New(n);
  if(!bnames) return NULL;
  for(i=0;i<n;i++){
      PyObject *val = Py_BuildValue("s",blocknames[i]);
      if(!val) { Py_DECREF(bnames); return NULL; }
      PyList_SET_ITEM(bnames,i,val);
  }
  free(blocknames);
  return bnames;
}

static PyObject *
pairix_get_chr1_col(PairixObject *self)
{
  return Py_BuildValue("i", ti_get_sc(self->tb->idx));
}

static PyObject *
pairix_get_chr2_col(PairixObject *self)
{
  return Py_BuildValue("i", ti_get_sc2(self->tb->idx));
}

static PyObject *
pairix_get_startpos1_col(PairixObject *self)
{
  return Py_BuildValue("i", ti_get_bc(self->tb->idx));
}

static PyObject *
pairix_get_startpos2_col(PairixObject *self)
{
  return Py_BuildValue("i", ti_get_bc2(self->tb->idx));
}

static PyObject *
pairix_get_endpos1_col(PairixObject *self)
{
  return Py_BuildValue("i", ti_get_ec(self->tb->idx));
}

static PyObject *
pairix_get_endpos2_col(PairixObject *self)
{
  return Py_BuildValue("i", ti_get_ec2(self->tb->idx));
}

static PyObject *
pairix_exists(PairixObject *self, PyObject *args)
{
    const char *key;
    if (!PyArg_ParseTuple(args, "s:exists", &key)){
        PyErr_SetString(PairixError, "Argument error! exists() takes the following args: <key_str>. Key_str is a str formatted as: '{CHR}' (1D) or  '{CHR}|{CHR}' (2D). (e.g. 'chr1|chr2')\n");
        return Py_BuildValue("i", -1);
    }
  return Py_BuildValue("i", ti_get_tid(self->tb->idx, key)!=-1?1:0); // returns 1 if key exists, 0 if not, -1 if usage is wrong.
}

static PyObject *
pairix_exists2(PairixObject *self, PyObject *args)
{
    char *chr1, *chr2, *sname;
    char region_split_character = get_region_split_character(self->tb);
    int h;
    if (!PyArg_ParseTuple(args, "ss:exists2", &chr1, &chr2)){
        PyErr_SetString(PairixError, "Argument error! exists2() takes the following args: <seqname1(chr1)> <seqname2(chr2)>\n");
        return Py_BuildValue("i", -1);
    }
    /* concatenate chromosomes */
    sname = (char*)malloc(strlen(chr1)+strlen(chr2)+2);
    strcpy(sname, chr1);
    h=strlen(sname);
    sname[h]= region_split_character;
    strcpy(sname+h+1, chr2);
    int res = ti_get_tid(self->tb->idx, sname)!=-1?1:0;
    free(sname);
    return Py_BuildValue("i", res); // returns 1 if key exists, 0 if not, -1 if usage is wrong.
}

static PyObject *
pairix_get_header(PairixObject *self)
{
    sequential_iter_t *siter;
    int len;
    int n=0; // number of header lines
    char *s;

    const ti_conf_t *pconf = ti_get_conf(self->tb->idx);
    siter = ti_query_general(self->tb, 0, 0, 0);
    while ((s = sequential_ti_read(siter, &len)) != 0) {
        if ((int)(*s) != pconf->meta_char) break;
        n++;
    }
    destroy_sequential_iter(siter);
    bgzf_seek(self->tb->fp, 0, SEEK_SET);

    PyObject *headers = PyList_New(n);
    if(!headers) return NULL;
    siter = ti_query_general(self->tb, 0, 0, 0);
    int i=0;
    while ((s = sequential_ti_read(siter, &len)) != 0) {
        if ((int)(*s) != pconf->meta_char) break;
        PyObject *val = Py_BuildValue("s", s);
        if(!val) { Py_DECREF(headers); return NULL; }
        PyList_SET_ITEM(headers, i, val);
        i++;
    }
    destroy_sequential_iter(siter);
    bgzf_seek(self->tb->fp, 0, SEEK_SET);
    return(headers);
}

static PyObject *
pairix_get_chromsize(PairixObject *self)
{
    sequential_iter_t *siter;
    int len;
    int n=0; // number of header lines corresponding to chromsize
    char *s;

    const ti_conf_t *pconf = ti_get_conf(self->tb->idx);
    siter = ti_query_general(self->tb, 0, 0, 0);
    while ((s = sequential_ti_read(siter, &len)) != 0) {
        if ((int)(*s) != pconf->meta_char) break;
        if(strncmp(s, "#chromsize: ", 12)==0)
            n++;
    }
    destroy_sequential_iter(siter);
    bgzf_seek(self->tb->fp, 0, SEEK_SET);

    PyObject *headers = PyList_New(n);
    if(!headers) return NULL;
    siter = ti_query_general(self->tb, 0, 0, 0);
    int i=0;
    while ((s = sequential_ti_read(siter, &len)) != 0) {
        if ((int)(*s) != pconf->meta_char) break;
        if(strncmp(s, "#chromsize: ", 12)==0) {
            PyObject *header_line = PyList_New(2);
            int j=12; do { j++; } while(s[j]!=' ' && s[j]!='\t');
            char b=s[j]; s[j]=0;
            PyObject *chr = Py_BuildValue("s", s + 12);
            if(!chr) { Py_DECREF(header_line); Py_DECREF(headers); return NULL; }
            PyList_SET_ITEM(header_line, 0, chr);
            s[j]=b;
            PyObject *val = Py_BuildValue("s", s + j + 1);
            if(!val) { Py_DECREF(header_line); Py_DECREF(headers);  return NULL; }
            PyList_SET_ITEM(header_line, 1, val);
            PyList_SET_ITEM(headers, i, header_line);
            i++;
        }
    }
    destroy_sequential_iter(siter);
    bgzf_seek(self->tb->fp, 0, SEEK_SET);
    return(headers);
}

static PyObject *
pairix_bgzf_block_count(PairixObject *self, PyObject *args)
{
    char *chr1, *chr2, *sname;
    char region_split_character = get_region_split_character(self->tb);
    int h;
    if (!PyArg_ParseTuple(args, "ss:bgzf_block_count", &chr1, &chr2)){
        PyErr_SetString(PairixError, "Argument error! bgzf_block_count() takes the following args: <seqname1(chr1)> <seqname2(chr2)>\n");
        return Py_BuildValue("i", -1);
    }
    /* concatenate chromosomes */
    sname = (char*)malloc(strlen(chr1)+strlen(chr2)+2);
    strcpy(sname, chr1);
    h=strlen(sname);
    sname[h]= region_split_character;
    strcpy(sname+h+1, chr2);
    int tid = ti_get_tid(self->tb->idx, sname);
    free(sname);
    if (tid == -1) {
        return Py_BuildValue("i", 0);
    }

    int res = get_nblocks(self->tb->idx, tid, self->tb->fp);
    return Py_BuildValue("i", res);
}


static PyObject *
pairix_check_triangle(PairixObject *self)
{
    int res = check_triangle(self->tb->idx);
    if(res>=0) {
        if(res == 1) printf("The file is a triangle.\n");
        else if(res == 0) printf("The file is not a triangle.\n");
        return Py_BuildValue("i", res);
    } else {
        if(res == -1) PyErr_SetString(PairixError, "Cannot retrieve seqnames.\n");
        else if(res == -2) PyErr_SetString(PairixError, "The file is 1D-indexed (option not applicable)\n");
        return Py_BuildValue("i", res);
    }
}

static PyObject *
pairix_repr(PairixObject *self)
{
#if PY_MAJOR_VERSION < 3
    return PyString_FromFormat("<pypairix fn=\"%s\">", self->fn);
#else
    return PyUnicode_FromFormat("<pypairix fn=\"%s\">", self->fn);
#endif
}

static PyMethodDef pairix_methods[] = {
    {
        "query",
        (PyCFunction)pairix_query,
        METH_VARARGS,
        PyDoc_STR("Retrieve items within a region.\n\n"
                  "    >>> tb.query(\"chr1\", 1000, 2000)\n"
                  "    <pypairix.iter at 0x17b86e50>\n\n"
                  "Parameters\n"
                  "----------\n"
                  "name : str\n"
                  "    Name of the sequence in the file.\n"
                  "start : int\n"
                  "    Start of the query region.\n"
                  "end : int\n"
                  "    End of the query region.\n")
    },
    {
        "queryi",
        (PyCFunction)pairix_queryi,
        METH_VARARGS,
        PyDoc_STR("Retrieve items within a region.\n\n"
                  "    >>> tb.queryi(0, 1000, 2000)\n"
                  "    <pypairix.iter at 0x17b86e50>\n\n"
                  "Parameters\n"
                  "----------\n"
                  "tid : int\n"
                  "    Index of the sequence in the file (first is 0).\n"
                  "start : int\n"
                  "    Start of the query region.\n"
                  "end : int\n"
                  "    End of the query region.\n")
    },
    {
        "querys",
        (PyCFunction)pairix_querys,
        METH_VARARGS,
        PyDoc_STR("Retrieve items within a region.\n\n"
                  "    >>> tb.querys(\"chr1:1000-2000\")\n"
                  "    <pypairix.iter at 0x17b86e50>\n\n"
                  "Parameters\n"
                  "----------\n"
                  "region : str\n"
                  "    Query string like \"seq:start-end\".\n")
    },
    {
        "query2D",
        (PyCFunction)pairix_query_2D,
        METH_VARARGS,
        PyDoc_STR("Retrieve items within a region.\n\n"
                  "    >>> pr.query2D(\"chr1\", 1000, 2000, \"chr2\", 1, 10000)\n"
                  "    <pypairix.iter at 0x17b86e50>\n\n"
                  "Parameters\n"
                  "----------\n"
                  "name : str\n"
                  "    Name of the sequence in the file.\n"
                  "start : int\n"
                  "    Start of the query region.\n"
                  "end : int\n"
                  "    End of the query region.\n"
                  "name2 : str\n"
                  "    Name of the sequence in the file.\n"
                  "start2 : int\n"
                  "    Start of the query region.\n"
                  "end2 : int\n"
                  "    End of the query region.\n"
                  "flip value: int\n"
                  "    If == 1, will attempt to flip chromsomes order on an error.\n")
    },
    {
        "queryi2D",
        (PyCFunction)pairix_queryi_2D,
        METH_VARARGS,
        PyDoc_STR("Retrieve items within a region.\n\n"
                  "    >>> pr.queryi(0, 1000, 2000)\n"
                  "    <pypairix.iter at 0x17b86e50>\n\n"
                  "Parameters\n"
                  "----------\n"
                  "tid : int\n"
                  "    Index of the sequence in the file (first is 0).\n"
                  "start : int\n"
                  "    Start of the query region.\n"
                  "end : int\n"
                  "    End of the query region.\n")
    },
    {
        "querys2D",
        (PyCFunction)pairix_querys_2D,
        METH_VARARGS,
        PyDoc_STR("Retrieve items within a region.\n\n"
                  "    >>> pr.querys2D(\"chr1:1000-2000|chr2:1-10000\")\n"
                  "    <pypairix.iter at 0x17b86e50>\n\n"
                  "Parameters\n"
                  "----------\n"
                  "region : str\n"
                  "    Query string like \"chr1:start-end|chr2:start-end\".\n"
                  "flip value: int\n"
                  "    If == 1, will attempt to flip chromsomes order on an error.\n")
    },
    {
       "get_linecount",
       (PyCFunction)pairix_get_linecount,
        METH_VARARGS,
        PyDoc_STR("Retrieve total linecount of the file.\n\n")
    },
    {
       "get_blocknames",
       (PyCFunction)pairix_get_blocknames,
        METH_VARARGS,
        PyDoc_STR("Retrieve list of keys (either chromosomes(1D-indexed) or chromosome pairs(2D-indexed)).\n\n")
    },
    {
       "get_chr1_col",
       (PyCFunction)pairix_get_chr1_col,
        METH_VARARGS,
        PyDoc_STR("Retrieve the 0-based column index of the first chromosome.\n\n")
    },
    {
       "get_chr2_col",
       (PyCFunction)pairix_get_chr2_col,
        METH_VARARGS,
        PyDoc_STR("Retrieve the 0-based column index of the second chromosome.\n\n")
    },
    {
       "get_startpos1_col",
       (PyCFunction)pairix_get_startpos1_col,
        METH_VARARGS,
        PyDoc_STR("Retrieve the 0-based column index of the start position on the first chromosome.\n\n")
    },
    {
       "get_startpos2_col",
       (PyCFunction)pairix_get_startpos2_col,
        METH_VARARGS,
        PyDoc_STR("Retrieve the 0-based column index of the start potision on the second chromosome.\n\n")
    },
    {
       "get_endpos1_col",
       (PyCFunction)pairix_get_endpos1_col,
        METH_VARARGS,
        PyDoc_STR("Retrieve the 0-based column index of the end position of the first chromosome.\n\n")
    },
    {
       "get_endpos2_col",
       (PyCFunction)pairix_get_endpos2_col,
        METH_VARARGS,
        PyDoc_STR("Retrieve the 0-based column index of the end position of the second chromosome.\n\n")
    },
    {
       "exists",
       (PyCFunction)pairix_exists,
        METH_VARARGS,
        PyDoc_STR("Check if key exists(1 if exists, 0 if not, -1 if wrong usage.)\n\n")
    },
    {
       "exists2",
       (PyCFunction)pairix_exists2,
        METH_VARARGS,
        PyDoc_STR("Check if key exists2(1 if exists, 0 if not, -1 if wrong usage.)\n\n")
    },
    {
       "get_header",
       (PyCFunction)pairix_get_header,
        METH_VARARGS,
        PyDoc_STR("return header as a list of strings\n\n")
    },
    {
       "get_chromsize",
       (PyCFunction)pairix_get_chromsize,
        METH_VARARGS,
        PyDoc_STR("return chromsize as a list of 2-element lists \
                   [chromosome_name, chromosome_size]\n\n")
    },
    {
       "bgzf_block_count",
       (PyCFunction)pairix_bgzf_block_count,
        METH_VARARGS,
        PyDoc_STR("return count of bgzf blocks for given key\n\n")
    },
    {
       "check_triangle",
       (PyCFunction)pairix_check_triangle,
        METH_VARARGS,
        PyDoc_STR("return 1 if triangle (a chromosome pair occurs only in one direction (e.g. if chr1|chr2 exists, chr2|chr1 doesn't), 0 if not, -1/-2 if error\n\n")
    },
    /*
    {
        "header",
        (PyCFunction)tabix_header,
        METH_VARARGS,
        PyDoc_STR("Get the header for a file (VCF, SAM, GTF, etc.).\n\n"
                  "    >>> tb.header()\n")
    },
    */
    {NULL, NULL}           /* sentinel */
};


static PyTypeObject Pairix_Type = {
    /* The ob_type field must be initialized in the module init function
     * to be portable to Windows without using C++. */
    PyVarObject_HEAD_INIT(NULL, 0)
    "pypairix.open",              /*tp_name*/
    sizeof(PairixObject),        /*tp_basicsize*/
    0,                          /*tp_itemsize*/
    /* methods */
    (destructor)pairix_dealloc,  /*tp_dealloc*/
    0,                          /*tp_print*/
    0,                          /*tp_getattr*/
    0,                          /*tp_setattr*/
    0,                          /*tp_compare*/
    (reprfunc)pairix_repr,       /*tp_repr*/
    0,                          /*tp_as_number*/
    0,                          /*tp_as_sequence*/
    0,                          /*tp_as_mapping*/
    0,                          /*tp_hash*/
    0,                          /*tp_call*/
    0,                          /*tp_str*/
    0,                          /*tp_getattro*/
    0,                          /*tp_setattro*/
    0,                          /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,         /*tp_flags*/
    0,                          /*tp_doc*/
    0,                          /*tp_traverse*/
    0,                          /*tp_clear*/
    0,                          /*tp_richcompare*/
    0,                          /*tp_weaklistoffset*/
    0,                          /*tp_iter*/
    0,                          /*tp_iternext*/
    pairix_methods,              /*tp_methods*/
    0,                          /*tp_members*/
    0,                          /*tp_getset*/
    0,                          /*tp_base*/
    0,                          /*tp_dict*/
    0,                          /*tp_descr_get*/
    0,                          /*tp_descr_set*/
    0,                          /*tp_dictoffset*/
    0,                          /*tp_init*/
    0,                          /*tp_alloc*/
    (newfunc)pairix_new,         /*tp_new*/
    0,                          /*tp_free*/
    0,                          /*tp_is_gc*/
};


static PyMethodDef pairix_functions[] = {
    {NULL, NULL} /* sentinel */
};

PyDoc_STRVAR(module_doc,
"Python interface to pairix.");

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef pypairixmodule = {
    PyModuleDef_HEAD_INIT,
    "pypairix",
    module_doc,
    -1,
    pairix_functions,
    NULL,
    NULL,
    NULL,
    NULL
};
#endif

#if PY_MAJOR_VERSION < 3
PyMODINIT_FUNC initpypairix(void)
#else
PyMODINIT_FUNC PyInit_pypairix(void)
#endif
{
    PyObject *m;

    if (PyType_Ready(&Pairix_Type) < 0)
        goto fail;
    if (PyType_Ready(&PairixIterator_Type) < 0)
        goto fail;

#if PY_MAJOR_VERSION < 3
    m = Py_InitModule3("pypairix", pairix_functions, module_doc);
#else
    m = PyModule_Create(&pypairixmodule);
#endif
    if (m == NULL)
        goto fail;

    if (PairixError == NULL) {
        PairixError = PyErr_NewException("pypairix.PairixError", NULL, NULL);
        if (PairixError == NULL)
            goto fail;
    }
    Py_INCREF(PairixError);
    PyModule_AddObject(m, "PairixError", PairixError);

    if (PairixWarning == NULL) {
        PairixWarning = PyErr_NewException("pypairix.PairixWarning", NULL, NULL);
        if (PairixWarning == NULL)
            goto fail;
    }
    Py_INCREF(PairixWarning);
    PyModule_AddObject(m, "PairixWarning", PairixWarning);

    PyModule_AddObject(m, "open", (PyObject *)&Pairix_Type);
    PyModule_AddObject(m, "iter", (PyObject *)&PairixIterator_Type);
    PyObject *module_dict, *indexer_func, *mod_name;
    mod_name = PYOBJECT_FROM_STRING("pypairix");
    module_dict = PyModule_GetDict(m);
    indexer_func = PyCFunction_NewEx(&indexer_methods[0], (PyObject*)NULL, mod_name);
    PyDict_SetItemString(module_dict, "build_index", indexer_func);
    PyDict_SetItemString(module_dict, "__version__", PYOBJECT_FROM_STRING(PACKAGE_VERSION));

#if PY_MAJOR_VERSION >= 3
    return m;
#endif

 fail:
#if PY_MAJOR_VERSION < 3
    return;
#else
    return NULL;
#endif
}
