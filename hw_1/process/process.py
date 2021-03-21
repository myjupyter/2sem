import os
import math

import pymorphy2

from collections import Counter

from .ru_stop_words import stop_words

morph = pymorphy2.MorphAnalyzer()

alphabet_upper = "А,Б,В,Г,Д,Е,Ё,Ж,З,И,Й,К,Л,М,Н,О,П,Р,С,Т,У,Ф,Х,Ц,Ч,Ш,Щ,Ъ,Ы,Ь,Э,Ю,Я"
alphabet_lower = alphabet_upper.lower()

alphabet = set(alphabet_upper.split(',')).union(alphabet_lower.split(','))

class Query:
    def __init__(self, query = None):
        if not isinstance(query, str):
            raise Exception('query should be a string')

        self.__query_original = query
        self.__query = clear_sentence(query)
        self.__counter = counter_from(self.__query)
   
    
    def tf(self, word):
        norm_word = normalize_word(word)
        return 0 if self.__counter.get(norm_word) is None else self.__counter[norm_word]
   
    @property
    def coords(self):
        return set(self.__counter.keys())

    def tf_vector(self, method=1):
        tf = {coord: self.tf(coord) for coord in self.coords}
        if method >= 1:
            return tf
        tf_max = max(tf.items(), key = lambda x: x[1])[1]
        if tf_max == 0:
            return {coord: 0.0 for coord in self.coords}
        return {coord: 0.4 + 0.6 * self.tf(coord) / tf_max for coord in self.coords}
        
class Documents:
    def __init__(self, doc_paths = None):
        if not isinstance(doc_paths, list):
            raise Exception('wrong type doc_path, doc_path should be a string')

        self.__sentences = []
        for doc_path in doc_paths:
            with open(doc_path) as file:
                text = file.read()

            sentences = text.split('.')
            name = os.path.basename(doc_path)
            sentences = [Sentence(sentence, name) for sentence in sentences]
            
            self.__sentences.extend(sentences)


    @property
    def counter(self):
        return self.__counter
   
    def idf(self, coords):
        return {coord: self.__idf(coord) for coord in coords}

    def __idf(self, word):
        norm_word = normalize_word(word)
        n_sentences = len(self.__sentences)
        
        has_sent = 0
        for sentence in self.__sentences:
            if sentence.counter.get(norm_word) is not None:
                has_sent += 1

        if has_sent == 0:
            has_sent += 1
        
        return math.log10(n_sentences / has_sent)
    
    def search_n(self, query, n, method=1):
        query = Query(query) 
        query_tf_vector = query.tf_vector(method)
        
        pairs = {}
        for i in range(len(self.__sentences)):
            tf = 0
            if method >= 1:
                tf = self.__sentences[i].tf()
            else:
                tf = self.__sentences[i].tf2()
            tfidf = tf_idf(tf, self.idf(query.coords))
            pairs[i] = tfidf

        
        cos_pairs = [] 
        for i, v in pairs.items():
            cos = cosine(query_tf_vector, v)
            print(self.__sentences[i], query_tf_vector, v, cos)
            cos_pairs.append((i, cos))
        
        
        l = n if len(cos_pairs) > n else len(cos_pairs) 
        cos_pairs = sorted(cos_pairs, key=lambda x: x[1], reverse=True)[:l]
        res = []
        for i, cos in cos_pairs:
            self.__sentences[i].set_weight(cos)
            res.append(self.__sentences[i])
        
        return res
    
            
    @property
    def counter(self):
        return self.__counter
    
class Sentence:
    def __init__(self, sentence, doc_name):
        self.__weight = 0 
        self.__doc_name = doc_name
        self.__original_sentence = sentence
        self.__sentence = clear_sentence(sentence)
        self.__counter = counter_from(self.__sentence) 

    def __str__(self):
        return self.__original_sentence

    def set_weight(self, weight):
        self.__weight = weight
    
    @property
    def weight(self):
        return self.__weight

    @property
    def doc_name(self):
        return self.__doc_name
    
    @property
    def counter(self):
        return self.__counter

    def tf(self):
        return self.__counter
        #return {coord: self.__tf(coord) for coord in coords}

    def tf2(self):
        tf = self.tf()
        if len(tf) == 0:
            return {}
        tf_max = max(tf.items(), key = lambda x: x[1])[1]
        if tf_max == 0:
            return {coord: 0.0 for coord in coords}
        return {coord: 0.4 + 0.6 * tf / tf_max for coord, tf in tf.items()}
    
    def __tf(self, word):
        norm_word = normalize_word(word)
        return 0 if self.__counter.get(norm_word) is None else self.__counter[norm_word]


def normalize_word(word):
    p = morph.parse(word)[0]
    return p.normal_form

def counter_from(sentence):
    c = Counter(sentence.split(' '))
    if c.get('') is not None:
        del c['']
    return c

def clear_sentence(sentence):
    sentence = sentence.lower()
    return ' '.join([ 
        check_stop_word(normalize_word(''.join([c if c in alphabet else ' ' for c in word]).strip()))
        for word in sentence.split(' ') 
    ]).strip()
    
def check_stop_word(word):
    if word in stop_words:
        return ' '
    return word

def cosine(tf_idf1, tf_idf2):
    norm = L2(tf_idf1) * L2(tf_idf2) 
    if norm == 0.0:
        return norm
    n = scal_mul(tf_idf1, tf_idf2) / norm 
#    print(tf_idf1, tf_idf2, n)
    return n
        
def L2(vector):
    return math.sqrt(sum([x*x for x in vector.values()]))  

def scal_mul(v1, v2):
    s = 0
    for coord in v1.keys():
        if v2.get(coord) is not None:
            s += v1[coord] * v2[coord]
    return s

def tf_idf(tf, idf):
    tfidf = {}
    for coord in tf.keys():
        if idf.get(coord) is not None:
            tfidf[coord] = tf[coord] * idf[coord]
        else:
            tfidf[coord] = 0
    return tfidf

def process(text):
    sentences = text.split('.')

    new_text = [clear_sentence(sentence) for sentence in sentences]
    return new_text 
