"""
================================================================================
   GENOMIC WORDS DISCOVERY TOOL v3.0
   Проверка гипотезы: осмысленные английские слова в генах человека
   Автор: Андрей Шимельфениг
   Дата: 2026
================================================================================

   Гипотеза: Гены, связанные с определенной функцией, содержат ДНК-последовательности,
             которые при 5-битном декодировании дают английские слова,
             семантически связанные с функцией гена.

   Пример: ASPA (миелин, скорость сигнала) → "LIGHT"
           SCN9A (натриевые каналы, заряд) → "ELECTRON"

   Статистика: p-value методом перестановок (10,000 итераций)
   Контроли: 4 типа контрольных экспериментов
================================================================================
"""

import requests
import json
import hashlib
import os
import random
import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
from collections import Counter
import time


# ==============================================================================
# ЧАСТЬ 1. 5-БИТНОЕ КОДИРОВАНИЕ (СЛОВА → ДНК)
# ==============================================================================

class Word2DNA:
    """
    5-bit encoding: A-Z + спецсимволы → ДНК
    Каждый символ = 5 бит → padding до 6 бит → 3 нуклеотида
    Длина ДНК = длина_слова × 3
    """

    # Таблица символов → 5-битный код
    CHAR_TO_BIN = {
        'A': '00000', 'B': '00001', 'C': '00010', 'D': '00011',
        'E': '00100', 'F': '00101', 'G': '00110', 'H': '00111',
        'I': '01000', 'J': '01001', 'K': '01010', 'L': '01011',
        'M': '01100', 'N': '01101', 'O': '01110', 'P': '01111',
        'Q': '10000', 'R': '10001', 'S': '10010', 'T': '10011',
        'U': '10100', 'V': '10101', 'W': '10110', 'X': '10111',
        'Y': '11000', 'Z': '11001',
        ' ': '11010', '.': '11011', ',': '11100',
        '-': '11101', '_': '11110', ':': '11111'
    }

    # Обратная таблица
    BIN_TO_CHAR = {v: k for k, v in CHAR_TO_BIN.items()}

    # 2 бита → нуклеотид
    PAIR_TO_NUC = {
        '00': 'A',  # Aденин
        '01': 'C',  # Цитозин
        '10': 'G',  # Гуанин
        '11': 'T'  # Тимин
    }

    # Нуклеотид → 2 бита
    NUC_TO_PAIR = {v: k for k, v in PAIR_TO_NUC.items()}

    @classmethod
    def word_to_dna(cls, word: str) -> str:
        """
        Преобразует английское слово в ДНК-последовательность.

        Алгоритм:
        1. Каждая буква → 5 бит
        2. Добавить padding (0) в конец для кратности 2
        3. Разбить на пары по 2 бита
        4. Каждая пара → нуклеотид (A,C,G,T)

        Длина: len(word) × 3 нуклеотида
        """
        word = word.upper()

        # 1. Конвертируем буквы в биты
        binary = ''
        for ch in word:
            if ch not in cls.CHAR_TO_BIN:
                raise ValueError(f"Недопустимый символ: '{ch}'. Используйте A-Z, пробел, ., , - _ :")
            binary += cls.CHAR_TO_BIN[ch]

        # 2. Padding для четности
        if len(binary) % 2 != 0:
            binary += '0'

        # 3. Пары → нуклеотиды
        dna = ''
        for i in range(0, len(binary), 2):
            pair = binary[i:i + 2]
            dna += cls.PAIR_TO_NUC.get(pair, 'N')

        return dna

    @classmethod
    def dna_to_word(cls, dna: str) -> str:
        """
        Декодирует ДНК обратно в слово (для проверки).
        """
        # 1. Нуклеотиды → бинарные пары
        binary = ''
        for nuc in dna:
            if nuc not in cls.NUC_TO_PAIR:
                return f"[?{nuc}?]"
            binary += cls.NUC_TO_PAIR[nuc]

        # 2. Удаляем padding (нули в конце)
        while binary.endswith('0'):
            binary = binary[:-1]

        # 3. Добиваем до кратности 5
        while len(binary) % 5 != 0:
            binary += '0'

        # 4. Разбиваем по 5 бит и декодируем
        word = ''
        for i in range(0, len(binary), 5):
            chunk = binary[i:i + 5]
            if chunk in cls.BIN_TO_CHAR:
                word += cls.BIN_TO_CHAR[chunk]
            else:
                word += '?'

        return word

    @classmethod
    def validate_encoding(cls, word: str) -> bool:
        """Проверяет, что кодирование/декодирование работает"""
        dna = cls.word_to_dna(word)
        decoded = cls.dna_to_word(dna)
        is_valid = decoded == word.upper()
        if not is_valid:
            print(f"❌ Ошибка: '{word}' → '{dna}' → '{decoded}'")
        return is_valid


# ==============================================================================
# ЧАСТЬ 2. ПОЛУЧЕНИЕ ДАННЫХ ИЗ UCSC
# ==============================================================================

class UCSCFetcher:
    """Получение геномных последовательностей с кэшированием"""

    def __init__(self, cache_dir: str = "./genome_cache"):
        self.cache_dir = cache_dir
        os.makedirs(cache_dir, exist_ok=True)
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'GenomicWords/3.0 (research; andrey.shimelfenig@example.com)'
        })

    def fetch_sequence(self, genome: str, chrom: str,
                       start: int, end: int,
                       use_cache: bool = True) -> Optional[str]:
        """
        Получает последовательность из UCSC.
        Координаты: 1-based, start включительно, end включительно.
        """
        # Ключ кэша
        cache_key = f"{genome}_{chrom}_{start}_{end}"
        cache_file = f"{self.cache_dir}/{hashlib.md5(cache_key.encode()).hexdigest()}.fasta"

        # Проверка кэша
        if use_cache and os.path.exists(cache_file):
            try:
                with open(cache_file, 'r') as f:
                    lines = f.readlines()
                    seq = ''.join(line.strip() for line in lines[1:])
                    return seq.upper()
            except:
                pass

        # Запрос к UCSC API
        url = f"https://api.genome.ucsc.edu/getData/sequence"
        params = {
            "genome": genome,
            "chrom": chrom,
            "start": start - 1,  # UCSC uses 0-based
            "end": end
        }

        try:
            resp = self.session.get(url, params=params, timeout=30)
            resp.raise_for_status()
            data = resp.json()

            # Извлекаем последовательность
            if 'seq' in data:
                seq = data['seq']
            elif 'dna' in data and isinstance(data['dna'], dict):
                seq = data['dna'].get('seq', '')
            else:
                raise ValueError(f"Неизвестный формат: {list(data.keys())}")

            # Сохраняем в кэш
            if use_cache and seq:
                with open(cache_file, 'w') as f:
                    f.write(f">{cache_key}\n")
                    for i in range(0, len(seq), 60):
                        f.write(seq[i:i + 60] + '\n')

            return seq.upper() if seq else None

        except Exception as e:
            print(f"  ⚠ Ошибка получения {chrom}:{start}-{end}: {e}")
            return None


# ==============================================================================
# ЧАСТЬ 3. БАЗА ЗНАНИЙ: ГЕНЫ И КЛЮЧЕВЫЕ СЛОВА
# ==============================================================================

@dataclass
class Gene:
    """Информация о гене-кандидате"""
    name: str
    chrom: str
    start: int
    end: int
    function: str
    keywords: List[str]  # Слова, которые МЫ ОЖИДАЕМ найти


# База данных генов и ожидаемых слов
GENES_DATABASE = [
    Gene(
        name="ASPA",
        chrom="chr10",
        start=100188400,
        end=100188418,
        function="Миелинизация, скорость проведения нервных сигналов",
        keywords=["LIGHT", "SPEED", "FAST", "SIGNAL", "NERVE", "MYELIN"]
    ),
    Gene(
        name="SCN9A",
        chrom="chr2",
        start=166210400,
        end=166210413,
        function="Натриевый канал, проведение электрического заряда",
        keywords=["ELECTRON", "CHARGE", "SODIUM", "CHANNEL", "CURRENT", "VOLTAGE"]
    ),
    Gene(
        name="COL1A1",
        chrom="chr17",
        start=50183200,
        end=50183212,
        function="Коллаген, плотность костей, опора, гравитация",
        keywords=["GRAVITY", "MASS", "WEIGHT", "BONE", "COLLAGEN", "SUPPORT"]
    ),
    Gene(
        name="MT-ND5",
        chrom="chrM",
        start=12337,
        end=12351,
        function="Митохондриальная энергия, квантовые эффекты",
        keywords=["PLANCK", "QUANTUM", "ENERGY", "MITOCHONDRIA", "POWER"]
    ),
    Gene(
        name="TRPV1",
        chrom="chr17",
        start=3642010,
        end=3642022,
        function="Терморецепция, температура, тепло",
        keywords=["HEAT", "HOT", "TEMPERATURE", "BOLTZMANN", "THERMAL"]
    ),
    Gene(
        name="KCNH2",
        chrom="chr7",
        start=150950100,
        end=150950112,
        function="Калиевые каналы сердца, тонкая структура",
        keywords=["FINE", "STRUCTURE", "HEART", "POTASSIUM", "CHANNEL"]
    ),
    Gene(
        name="HBA1",
        chrom="chr16",
        start=176680,
        end=176692,
        function="Гемоглобин, перенос кислорода, количество частиц",
        keywords=["AVOGADRO", "PARTICLE", "COUNT", "NUMBER", "OXYGEN", "BLOOD"]
    ),
]


# ==============================================================================
# ЧАСТЬ 4. ПОИСК СОВПАДЕНИЙ И СТАТИСТИКА
# ==============================================================================

class GenomicWordSearcher:
    """Поиск слов в геноме со статистической оценкой"""

    def __init__(self, fetcher: UCSCFetcher):
        self.fetcher = fetcher
        self.results = []
        self.encoder = Word2DNA

        # Проверяем кодирование при инициализации
        self._validate_encoder()

    def _validate_encoder(self):
        """Проверяет, что кодирование работает правильно"""
        test_words = ["LIGHT", "DNA", "TEST", "HELLO"]
        for word in test_words:
            if not self.encoder.validate_encoding(word):
                raise RuntimeError(f"Кодировщик не прошел проверку на слове '{word}'")
        print("✅ Кодировщик прошел проверку")

    def hamming_similarity(self, seq1: str, seq2: str) -> float:
        """Нормализованное расстояние Хэмминга"""
        if len(seq1) != len(seq2):
            return 0.0
        matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
        return matches / len(seq1)

    def sliding_window_search(self, sequence: str, target: str) -> Tuple[float, int, str]:
        """Поиск наилучшего совпадения методом скользящего окна"""
        window_len = len(target)
        best_score = 0.0
        best_pos = -1
        best_match = ''

        for i in range(len(sequence) - window_len + 1):
            window = sequence[i:i + window_len]
            score = self.hamming_similarity(window, target)

            if score > best_score:
                best_score = score
                best_pos = i
                best_match = window

                if score == 1.0:  # Идеальное совпадение
                    break

        return best_score, best_pos, best_match

    def permutation_test(self, target: str, observed_score: float,
                         n_permutations: int = 10000) -> Dict:
        """
        Перестановочный тест для оценки значимости.

        Возвращает:
        - p_value: статистическая значимость
        - mean_random: среднее случайных совпадений
        - std_random: стандартное отклонение
        - z_score: Z-оценка
        """
        # Генерируем случайные последовательности с тем же GC-составом
        gc_content = sum(1 for nuc in target if nuc in 'GC') / len(target)
        nucleotides = ['A', 'C', 'G', 'T']
        probs = [
            (1 - gc_content) / 2,  # A
            gc_content / 2,  # C
            gc_content / 2,  # G
            (1 - gc_content) / 2  # T
        ]

        random_scores = []
        for _ in range(n_permutations):
            random_seq = ''.join(np.random.choice(nucleotides, size=len(target), p=probs))
            score = self.hamming_similarity(random_seq, target)
            random_scores.append(score)

        # Расчет статистик
        mean_random = np.mean(random_scores)
        std_random = np.std(random_scores)
        p_value = sum(1 for s in random_scores if s >= observed_score) / n_permutations
        z_score = (observed_score - mean_random) / std_random if std_random > 0 else 0

        return {
            'p_value': p_value,
            'mean_random': mean_random,
            'std_random': std_random,
            'z_score': z_score,
            'n_permutations': n_permutations
        }

    def search_word_in_gene(self, word: str, gene: Gene,
                            genome: str = 'hg38') -> Optional[Dict]:
        """
        Ищет слово в гене и оценивает статистическую значимость.
        """
        print(f"\n{'=' * 70}")
        print(f"🔍 ПОИСК: '{word}' в гене {gene.name}")
        print(f"   {gene.function}")
        print(f"{'=' * 70}")

        # 1. Кодируем слово в ДНК
        try:
            target_dna = self.encoder.word_to_dna(word)
        except ValueError as e:
            print(f"❌ Ошибка кодирования: {e}")
            return None

        print(f"\n📝 Слово: '{word}'")
        print(f"🧬 ДНК:   {target_dna}")
        print(f"📏 Длина: {len(target_dna)} нуклеотидов")
        print(f"🔄 Проверка: '{word}' → '{self.encoder.dna_to_word(target_dna)}'")

        # 2. Получаем последовательность гена
        print(f"\n📡 Запрос к UCSC: {gene.chrom}:{gene.start}-{gene.end}")
        gene_seq = self.fetcher.fetch_sequence(
            genome=genome,
            chrom=gene.chrom,
            start=gene.start,
            end=gene.end
        )

        if not gene_seq:
            print(f"❌ Не удалось получить последовательность гена")
            return None

        print(f"✅ Получено: {len(gene_seq)} нуклеотидов")
        print(f"   {gene_seq[:50]}...")

        # 3. Поиск совпадений
        print(f"\n🔎 Поиск совпадений...")
        best_score, best_pos, best_match = self.sliding_window_search(gene_seq, target_dna)

        # 4. Статистическая оценка
        print(f"\n📊 Статистический анализ...")
        stats = self.permutation_test(target_dna, best_score, n_permutations=10000)

        # 5. Формируем результат
        result = {
            'word': word,
            'gene': gene.name,
            'chrom': gene.chrom,
            'position': gene.start + best_pos,
            'target_dna': target_dna,
            'best_match': best_match,
            'similarity': best_score,
            'matches': int(best_score * len(target_dna)),
            'length': len(target_dna),
            'p_value': stats['p_value'],
            'z_score': stats['z_score'],
            'significant_05': stats['p_value'] < 0.05,
            'significant_01': stats['p_value'] < 0.01,
            'significant_001': stats['p_value'] < 0.001,
            'significant_0001': stats['p_value'] < 0.0001,
            'timestamp': time.time()
        }

        # 6. Вывод результатов
        print(f"\n{'=' * 70}")
        print(f"🎯 РЕЗУЛЬТАТ:")
        print(f"{'=' * 70}")
        print(f"   Слово:          {result['word']}")
        print(f"   Ген:            {result['gene']}")
        print(f"   Позиция:        {result['chrom']}:{result['position']}")
        print(f"   Целевая ДНК:    {result['target_dna']}")
        print(f"   Найденo:        {result['best_match']}")
        print(f"   Совпадений:     {result['matches']}/{result['length']} ({result['similarity']:.1%})")
        print(f"   P-value:        {result['p_value']:.6f}")
        print(f"   Z-score:        {result['z_score']:.2f}")

        # Звездочки значимости
        if result['significant_0001']:
            print(f"   ⭐⭐⭐⭐ p < 0.0001 - ЧРЕЗВЫЧАЙНО ЗНАЧИМО!")
        elif result['significant_001']:
            print(f"   ⭐⭐⭐ p < 0.001 - ОЧЕНЬ ЗНАЧИМО!")
        elif result['significant_01']:
            print(f"   ⭐⭐ p < 0.01 - ЗНАЧИМО!")
        elif result['significant_05']:
            print(f"   ⭐ p < 0.05 - НА ГРАНИ ЗНАЧИМОСТИ")
        else:
            print(f"   ❌ НЕ ЗНАЧИМО")

        self.results.append(result)
        return result

    def analyze_gene_hypothesis(self, gene: Gene, genome: str = 'hg38') -> List[Dict]:
        """
        Проверяет все гипотезы для данного гена.
        Ищет все ключевые слова в гене.
        """
        print(f"\n{'#' * 80}")
        print(f"# АНАЛИЗ ГЕНА: {gene.name}")
        print(f"# {gene.function}")
        print(f"# Ключевые слова: {', '.join(gene.keywords)}")
        print(f"{'#' * 80}")

        results = []
        for word in gene.keywords:
            result = self.search_word_in_gene(word, gene, genome)
            if result:
                results.append(result)

        return results

    def run_full_analysis(self, genome: str = 'hg38') -> pd.DataFrame:
        """
        Запускает полный анализ всех генов и всех ключевых слов.
        """
        print("\n" + "=" * 80)
        print("🧬 ПОЛНЫЙ ГЕНОМНЫЙ АНАЛИЗ")
        print("=" * 80)
        print(f"Генов для анализа: {len(GENES_DATABASE)}")
        print(f"Всего гипотез: {sum(len(g.keywords) for g in GENES_DATABASE)}")
        print("=" * 80)

        all_results = []
        for gene in GENES_DATABASE:
            results = self.analyze_gene_hypothesis(gene, genome)
            all_results.extend(results)

        # Сортируем по значимости
        all_results.sort(key=lambda x: x['p_value'])

        return pd.DataFrame(all_results)

    def print_summary_report(self, df: pd.DataFrame):
        """Печатает итоговый отчет"""
        print("\n" + "=" * 80)
        print("📊 ИТОГОВЫЙ ОТЧЕТ")
        print("=" * 80)

        if len(df) == 0:
            print("❌ Нет результатов")
            return

        # Значимые результаты
        significant = df[df['significant_01']].copy()

        print(f"\n📈 Всего проверено гипотез: {len(df)}")
        print(f"✅ Статистически значимых (p < 0.01): {len(significant)}")
        print(f"   Из них p < 0.001: {len(df[df['significant_001']])}")
        print(f"        p < 0.0001: {len(df[df['significant_0001']])}")

        if len(significant) > 0:
            print(f"\n{'⭐' * 40}")
            print(f"🌟 ЗНАЧИМЫЕ ОТКРЫТИЯ:")
            print(f"{'⭐' * 40}")

            for _, row in significant.iterrows():
                stars = ""
                if row['p_value'] < 0.0001:
                    stars = "⭐⭐⭐⭐"
                elif row['p_value'] < 0.001:
                    stars = "⭐⭐⭐"
                elif row['p_value'] < 0.01:
                    stars = "⭐⭐"

                print(f"\n  {stars} {row['word']} → {row['gene']}")
                print(f"     Совпадение: {row['matches']}/{row['length']} ({row['similarity']:.1%})")
                print(f"     Позиция: {row['chrom']}:{row['position']}")
                print(f"     p = {row['p_value']:.6f}, Z = {row['z_score']:.2f}")
                print(f"     ДНК: {row['target_dna']}")
                print(f"     Найдено: {row['best_match']}")

        else:
            print(f"\n❌ Значимых результатов не обнаружено.")
            print(f"   Максимальная точность: {df['similarity'].max():.1%}")
            print(f"   Минимальный p-value: {df['p_value'].min():.6f}")


# ==============================================================================
# ЧАСТЬ 5. КОНТРОЛЬНЫЕ ЭКСПЕРИМЕНТЫ
# ==============================================================================

class ControlExperiments:
    """Набор контрольных экспериментов для проверки специфичности"""

    def __init__(self, searcher: GenomicWordSearcher):
        self.searcher = searcher

    def control_1_random_words(self, gene: Gene, n_random: int = 100):
        """
        Контроль 1: Случайные слова той же длины.

        Ожидание: реальное слово (LIGHT, ELECTRON и т.д.) должно давать
        значительно более высокий score, чем случайные слова.
        """
        print(f"\n{'=' * 70}")
        print(f"🔬 КОНТРОЛЬ 1: Случайные слова vs реальное слово")
        print(f"{'=' * 70}")

        # Берем реальное слово из гипотезы
        real_word = gene.keywords[0]
        real_result = self.searcher.search_word_in_gene(real_word, gene)
        if not real_result:
            return

        real_score = real_result['similarity']

        # Генерируем случайные слова
        alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        random_scores = []

        for i in range(n_random):
            word_len = len(real_word)
            random_word = ''.join(random.choice(alphabet) for _ in range(word_len))
            result = self.searcher.search_word_in_gene(random_word, gene)
            if result:
                random_scores.append(result['similarity'])

        # Статистика
        mean_random = np.mean(random_scores)
        std_random = np.std(random_scores)
        percentile = sum(1 for s in random_scores if s >= real_score) / n_random * 100

        print(f"\n📊 РЕЗУЛЬТАТ КОНТРОЛЯ:")
        print(f"   Реальное слово '{real_word}': {real_score:.1%}")
        print(f"   Случайные слова: среднее {mean_random:.1%} ± {std_random:.1%}")
        print(f"   Процентиль реального слова: {100 - percentile:.1f}%")
        print(f"   Вывод: {'✅ СПЕЦИФИЧНО' if real_score > mean_random + 2 * std_random else '❌ НЕ СПЕЦИФИЧНО'}")

    def control_2_different_genes(self, word: str, correct_gene: Gene,
                                  other_genes: List[Gene]):
        """
        Контроль 2: То же слово в других генах.

        Ожидание: слово должно специфично находиться в "своем" гене,
        а не во всех подряд.
        """
        print(f"\n{'=' * 70}")
        print(f"🔬 КОНТРОЛЬ 2: Специфичность к гену")
        print(f"{'=' * 70}")
        print(f"   Слово: '{word}'")
        print(f"   Правильный ген: {correct_gene.name}")

        # Ищем в правильном гене
        correct_result = self.searcher.search_word_in_gene(word, correct_gene)
        if not correct_result:
            return

        correct_score = correct_result['similarity']

        # Ищем в других генах
        other_scores = []
        for gene in other_genes:
            if gene.name != correct_gene.name:
                result = self.searcher.search_word_in_gene(word, gene)
                if result:
                    other_scores.append(result['similarity'])

        print(f"\n📊 РЕЗУЛЬТАТ КОНТРОЛЯ:")
        print(f"   {correct_gene.name}: {correct_score:.1%}")
        for i, score in enumerate(other_scores[:3]):
            print(f"   Другой ген {i + 1}: {score:.1%}")

        is_specific = all(score < correct_score for score in other_scores)
        print(f"   Вывод: {'✅ СПЕЦИФИЧНО' if is_specific else '❌ НЕ СПЕЦИФИЧНО'}")

    def control_3_permuted_sequence(self, word: str, gene: Gene, n_permutations: int = 100):
        """
        Контроль 3: Перестановка нуклеотидов.

        Ожидание: сигнал зависит от порядка, а не от состава.
        """
        print(f"\n{'=' * 70}")
        print(f"🔬 КОНТРОЛЬ 3: Перестановка последовательности")
        print(f"{'=' * 70}")

        # Оригинальный поиск
        original_result = self.searcher.search_word_in_gene(word, gene)
        if not original_result:
            return

        original_score = original_result['similarity']
        target_dna = original_result['target_dna']

        # Получаем последовательность гена
        gene_seq = self.searcher.fetcher.fetch_sequence('hg38', gene.chrom, gene.start, gene.end)

        # Переставляем нуклеотиды
        permuted_scores = []
        for _ in range(n_permutations):
            seq_list = list(gene_seq)
            random.shuffle(seq_list)
            permuted_seq = ''.join(seq_list)

            score, _, _ = self.searcher.sliding_window_search(permuted_seq, target_dna)
            permuted_scores.append(score)

        mean_permuted = np.mean(permuted_scores)
        std_permuted = np.std(permuted_scores)

        print(f"\n📊 РЕЗУЛЬТАТ КОНТРОЛЯ:")
        print(f"   Оригинальная последовательность: {original_score:.1%}")
        print(f"   Переставленная: {mean_permuted:.1%} ± {std_permuted:.1%}")
        print(f"   Разница: {original_score - mean_permuted:.1%}")
        print(
            f"   Вывод: {'✅ ЗАВИСИТ ОТ ПОРЯДКА' if original_score > mean_permuted + 2 * std_permuted else '❌ ТОЛЬКО СОСТАВ'}")

    def run_all_controls(self):
        """Запускает все контрольные эксперименты"""
        print("\n" + "=" * 80)
        print("🧪 КОНТРОЛЬНЫЕ ЭКСПЕРИМЕНТЫ")
        print("=" * 80)

        # Контроль 1: ASPA и слово LIGHT
        aspa = next(g for g in GENES_DATABASE if g.name == "ASPA")
        self.control_1_random_words(aspa)

        # Контроль 2: LIGHT в других генах
        self.control_2_different_genes("LIGHT", aspa, GENES_DATABASE)

        # Контроль 3: Перестановка
        self.control_3_permuted_sequence("LIGHT", aspa)

        print("\n" + "=" * 80)
        print("✅ КОНТРОЛЬНЫЕ ЭКСПЕРИМЕНТЫ ЗАВЕРШЕНЫ")
        print("=" * 80)


# ==============================================================================
# ЧАСТЬ 6. ОСНОВНАЯ ПРОГРАММА
# ==============================================================================

def main():
    """Главная функция запуска"""

    print("""
    ╔══════════════════════════════════════════════════════════════╗
    ║     🧬 GENOMIC WORDS DISCOVERY TOOL v3.0                     ║
    ║     Проверка гипотезы: осмысленные слова в генах человека   ║
    ║     Автор: Андрей Шимельфениг                               ║
    ║     Дата: 2026                                              ║
    ╚══════════════════════════════════════════════════════════════╝
    """)

    # 1. Инициализация
    print("\n[1/5] Инициализация компонентов...")
    fetcher = UCSCFetcher(cache_dir="./genome_cache")
    searcher = GenomicWordSearcher(fetcher)

    # 2. Проверка кодирования
    print("\n[2/5] Проверка 5-битного кодировщика...")
    test_words = ["LIGHT", "DNA", "HELLO", "WORLD"]
    for word in test_words:
        dna = Word2DNA.word_to_dna(word)
        decoded = Word2DNA.dna_to_word(dna)
        print(f"   {word:8} → {dna:20} → {decoded}")

    # 3. Полный анализ
    print("\n[3/5] Запуск полного геномного анализа...")
    import pandas as pd
    results_df = searcher.run_full_analysis()

    # 4. Итоговый отчет
    print("\n[4/5] Формирование отчета...")
    searcher.print_summary_report(results_df)

    # 5. Контрольные эксперименты
    print("\n[5/5] Запуск контрольных экспериментов...")
    controls = ControlExperiments(searcher)
    controls.run_all_controls()

    # Сохраняем результаты
    if len(results_df) > 0:
        filename = f"genomic_words_results_{time.strftime('%Y%m%d_%H%M%S')}.csv"
        results_df.to_csv(filename, index=False)
        print(f"\n💾 Результаты сохранены в: {filename}")

    print("\n" + "=" * 80)
    print("🎉 АНАЛИЗ ЗАВЕРШЕН!")
    print("=" * 80)

    return results_df


if __name__ == "__main__":
    try:
        import pandas as pd
    except ImportError:
        print("📦 Устанавливаем pandas...")
        import subprocess

        subprocess.check_call(['pip', 'install', 'pandas'])
        import pandas as pd

    results = main()