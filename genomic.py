"""
Genomic Constants Mapping Tool
Author: Andrey Shimelfenig
Date: 2026
Version: 2.1 (Без зависимостей от scipy)
"""

import requests
import json
import numpy as np
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from collections import defaultdict
import hashlib
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import os
import random

@dataclass
class PhysicalConstant:
    """Физическая константа с метаданными"""
    name: str
    value: str
    gene: str
    chrom: str
    start: int
    end: int
    original_sequence: str

    @property
    def length(self) -> int:
        return len(self.original_sequence)

class DNACodingScheme:
    """
    Схема кодирования цифр в нуклеотиды.
    """

    # Основная схема: бинарное кодирование
    BINARY_MAP = {
        '0': '0000', '1': '0001', '2': '0010', '3': '0011',
        '4': '0100', '5': '0101', '6': '0110', '7': '0111',
        '8': '1000', '9': '1001', '.': '1010', '-': '1011',
        '+': '1100', 'e': '1101'
    }

    # Схема трансляции бинарных пар в нуклеотиды
    PAIR_TO_NUC = {
        '00': 'A',
        '01': 'C',
        '10': 'G',
        '11': 'T'
    }

    @classmethod
    def digit_to_sequence(cls, digit_str: str) -> str:
        """
        Преобразует строку с цифрами в последовательность нуклеотидов.

        Args:
            digit_str: Строка с цифрами (например, "299792458")

        Returns:
            Последовательность нуклеотидов
        """
        # 1. Конвертируем каждый символ в 4-битный код
        binary_str = ''
        for ch in digit_str:
            binary_str += cls.BINARY_MAP.get(ch, '0000')

        # 2. Разбиваем на пары
        pairs = [binary_str[i:i+2] for i in range(0, len(binary_str), 2)]

        # 3. Конвертируем пары в нуклеотиды
        sequence = ''
        for pair in pairs:
            sequence += cls.PAIR_TO_NUC.get(pair, 'N')

        return sequence

    @classmethod
    def generate_random_scheme(cls, seed: int = 42) -> Dict:
        """Генерирует случайную схему кодирования для контроля"""
        random.seed(seed)
        digits = ['0','1','2','3','4','5','6','7','8','9','.',',','-','+','e']
        pairs = ['00','01','02','03','04','05','06','07','08','09',
                '10','11','12','13','14','15']
        random.shuffle(pairs)
        return {d: pairs[i] for i, d in enumerate(digits)}

class GenomicSequenceFetcher:
    """Получение последовательностей из UCSC с кэшированием"""

    def __init__(self, cache_dir: str = "./cache"):
        self.cache_dir = cache_dir
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'GenomicConstants/2.0 (research; andrey.shimelfenig@example.com)'
        })
        # Создаём директорию кэша, если её нет
        os.makedirs(self.cache_dir, exist_ok=True)

    def fetch_sequence(self, genome: str, chrom: str,
                      start: int, end: int,
                      use_cache: bool = True) -> Optional[str]:
        """
        Получает последовательность с кэшированием.
        Координаты: 1-based, start включительно, end включительно.
        """
        cache_key = f"{genome}_{chrom}_{start}_{end}"
        cache_file = f"{self.cache_dir}/{hashlib.md5(cache_key.encode()).hexdigest()}.fasta"

        # Проверка кэша
        if use_cache and os.path.exists(cache_file):
            try:
                with open(cache_file, 'r') as f:
                    lines = f.readlines()
                    # Пропускаем первую строку (заголовок)
                    seq = ''.join(line.strip() for line in lines[1:])
                    print(f"  ✓ Загружено из кэша: {cache_key}")
                    return seq.upper()
            except Exception as e:
                print(f"  ⚠ Ошибка чтения кэша: {e}")

        # UCSC API (0-based, start включительно, end исключено)
        url = (f"https://api.genome.ucsc.edu/getData/sequence"
               f"?genome={genome}&chrom={chrom}"
               f"&start={start-1}&end={end}")

        print(f"  ↻ Запрос: {chrom}:{start}-{end}")

        try:
            resp = self.session.get(url, timeout=30)
            resp.raise_for_status()
            data = resp.json()

            # Извлечение последовательности из разных форматов ответа
            if 'dna' in data and isinstance(data['dna'], dict):
                seq = data['dna'].get('seq', '')
            elif 'seq' in data:
                seq = data['seq']
            else:
                raise ValueError(f"Неизвестный формат ответа: {list(data.keys())}")

            # Сохранение в кэш
            if use_cache and seq:
                with open(cache_file, 'w') as f:
                    f.write(f">{cache_key}\n")
                    for i in range(0, len(seq), 60):
                        f.write(seq[i:i+60] + "\n")
                print(f"  ✓ Сохранено в кэш: {cache_file}")

            return seq.upper() if seq else None

        except requests.exceptions.Timeout:
            print(f"  ✗ Таймаут: {chrom}")
            return None
        except requests.exceptions.HTTPError as e:
            print(f"  ✗ HTTP ошибка {e.response.status_code}")
            return None
        except Exception as e:
            print(f"  ✗ Ошибка: {e}")
            return None

class ConstantsMatcher:
    """
    Поиск соответствий констант в геноме со статистической обработкой.
    """

    def __init__(self, fetcher: GenomicSequenceFetcher):
        self.fetcher = fetcher
        self.results = []

    def calculate_similarity(self, seq1: str, seq2: str) -> float:
        """Расчёт расстояния Хэмминга (нормализованное)"""
        if len(seq1) != len(seq2):
            raise ValueError(f"Длины последовательностей не совпадают: {len(seq1)} vs {len(seq2)}")

        matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
        return matches / len(seq1)

    def generate_random_sequences(self, length: int,
                                 n: int = 1000,
                                 gc_content: float = 0.41) -> List[str]:
        """
        Генерация случайных последовательностей с заданным GC-составом.
        """
        sequences = []
        nucs = ['A', 'T', 'G', 'C']
        # Вероятности для GC-состава 41% (G+C=0.41, A+T=0.59)
        p_g = gc_content / 2
        p_c = gc_content / 2
        p_a = (1 - gc_content) / 2
        p_t = (1 - gc_content) / 2
        probs = [p_a, p_t, p_g, p_c]

        for _ in range(n):
            seq = ''.join(np.random.choice(nucs, size=length, p=probs))
            sequences.append(seq)

        return sequences

    def calculate_p_value(self, target_seq: str, observed_score: float,
                         n_permutations: int = 1000) -> Tuple[float, float, float]:
        """
        Вычисление p-value методом перестановок.

        Returns:
            (p_value, mean_random, std_random)
        """
        # Генерация случайных последовательностей
        random_seqs = self.generate_random_sequences(
            length=len(target_seq),
            n=n_permutations
        )

        # Расчёт случайных оценок
        random_scores = []
        for rand_seq in random_seqs:
            score = self.calculate_similarity(rand_seq, target_seq)
            random_scores.append(score)

        # Эмпирическое p-value
        n_extreme = sum(1 for score in random_scores if score >= observed_score)
        p_value = n_extreme / n_permutations

        # Среднее и стандартное отклонение
        mean_random = np.mean(random_scores)
        std_random = np.std(random_scores)

        return p_value, mean_random, std_random

    def sliding_window_search(self, sequence: str, target: str) -> Tuple[float, int, str]:
        """
        Поиск наилучшего совпадения методом скользящего окна.

        Returns:
            (best_score, best_position, best_match)
        """
        window_len = len(target)
        best_score = 0.0
        best_position = -1
        best_match = ''

        for i in range(len(sequence) - window_len + 1):
            window = sequence[i:i+window_len]
            score = self.calculate_similarity(window, target)

            if score > best_score:
                best_score = score
                best_position = i
                best_match = window

                # Если нашли идеальное совпадение, дальше можно не искать
                if score == 1.0:
                    break

        return best_score, best_position, best_match

    def analyze_constant(self, constant: PhysicalConstant,
                        debug: bool = False) -> Dict:
        """
        Полный анализ одной константы.
        """
        print(f"\n🔬 Анализ: {constant.name} ({constant.value})")

        # 1. Получаем последовательность гена
        print(f"  📥 Получение последовательности {constant.chrom}:{constant.start}-{constant.end}")
        gene_seq = self.fetcher.fetch_sequence(
            genome='hg38',
            chrom=constant.chrom,
            start=constant.start,
            end=constant.end
        )

        if not gene_seq:
            print(f"  ✗ Ошибка: не удалось получить последовательность")
            return {'error': 'Failed to fetch sequence', 'constant': constant.name}

        print(f"  ✓ Длина гена: {len(gene_seq)} п.н.")

        # 2. Кодируем константу в нуклеотиды
        target_seq = DNACodingScheme.digit_to_sequence(constant.value)
        print(f"  🔢 Константа '{constant.value}' → {target_seq}")
        print(f"  📏 Длина паттерна: {len(target_seq)} нуклеотидов")

        # 3. Поиск в гене
        best_score, best_position, best_match = self.sliding_window_search(
            gene_seq, target_seq
        )

        print(f"  🎯 Лучшее совпадение: {best_match}")
        print(f"  📍 Позиция: {constant.start + best_position}")
        print(f"  📊 Точность: {best_score:.3f} ({best_score*100:.1f}%)")

        # 4. Статистическая значимость
        print(f"  🧮 Расчёт p-value (1000 перестановок)...")
        p_value, mean_random, std_random = self.calculate_p_value(
            target_seq, best_score, n_permutations=1000
        )

        # Z-оценка
        z_score = (best_score - mean_random) / std_random if std_random > 0 else 0

        # Оценка значимости
        if p_value < 0.001:
            sig_mark = "***"
        elif p_value < 0.01:
            sig_mark = "**"
        elif p_value < 0.05:
            sig_mark = "*"
        else:
            sig_mark = "н.з."

        print(f"  📉 Случайные совпадения: μ={mean_random:.3f}, σ={std_random:.3f}")
        print(f"  📈 P-value: {p_value:.4f} {sig_mark}")
        print(f"  📊 Z-score: {z_score:.2f}")

        result = {
            'constant': constant.name,
            'value': constant.value,
            'target_sequence': target_seq,
            'best_match': best_match,
            'position': constant.start + best_position,
            'similarity': best_score,
            'p_value': p_value,
            'z_score': z_score,
            'significant': p_value < 0.05,
            'gene': constant.gene,
            'chrom': constant.chrom,
            'gene_region': f"{constant.start}-{constant.end}"
        }

        self.results.append(result)
        return result

    def generate_report(self, output_prefix: str = 'analysis_results'):
        """Генерация полного отчёта"""
        print(f"\n{'='*70}")
        print(f"📊 ИТОГОВЫЙ ОТЧЁТ")
        print(f"{'='*70}")

        if not self.results:
            print("❌ Нет результатов для отображения")
            return pd.DataFrame()

        # Создаём DataFrame
        df = pd.DataFrame(self.results)

        # Фильтрация значимых результатов
        significant = df[df['p_value'] < 0.05].copy()
        highly_significant = df[df['p_value'] < 0.01].copy()

        print(f"\n📈 Всего проанализировано констант: {len(df)}")
        print(f"✅ Значимых (p < 0.05): {len(significant)}")
        print(f"🔬 Высоко значимых (p < 0.01): {len(highly_significant)}")

        if len(significant) > 0:
            print(f"\n🌟 ЗНАЧИМЫЕ НАХОДКИ:")
            for _, row in significant.iterrows():
                stars = "***" if row['p_value'] < 0.001 else "**" if row['p_value'] < 0.01 else "*"
                print(f"\n  {row['constant']}:")
                print(f"    Точность: {row['similarity']:.3f} ({row['similarity']*100:.1f}%)")
                print(f"    P-value: {row['p_value']:.4f} {stars}")
                print(f"    Позиция: {row['chrom']}:{row['position']}")
                print(f"    Паттерн: {row['target_sequence']}")
                print(f"    Совпадение: {row['best_match']}")

        # Сохранение результатов
        json_file = f"{output_prefix}.json"
        csv_file = f"{output_prefix}.csv"

        df.to_json(json_file, orient='records', indent=2, force_ascii=False)
        df.to_csv(csv_file, index=False, encoding='utf-8')

        print(f"\n💾 Результаты сохранены:")
        print(f"  - {json_file}")
        print(f"  - {csv_file}")

        return df

    def plot_results(self, df: pd.DataFrame, output_file: str = 'analysis_plot.png'):
        """Визуализация результатов"""
        try:
            plt.figure(figsize=(12, 6))

            # График 1: Точность совпадений
            plt.subplot(1, 2, 1)
            constants = df['constant'].tolist()
            scores = df['similarity'].tolist()
            colors = ['red' if p < 0.05 else 'gray' for p in df['p_value']]

            bars = plt.bar(range(len(constants)), scores, color=colors, alpha=0.7)
            plt.axhline(y=0.5, color='black', linestyle='--', alpha=0.5, label='Случайный уровень')
            plt.xticks(range(len(constants)), constants, rotation=45, ha='right')
            plt.ylabel('Точность совпадения')
            plt.title('Точность совпадений с константами')
            plt.legend()

            # Добавим значения на столбцы
            for i, (bar, score) in enumerate(zip(bars, scores)):
                plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                        f'{score:.2f}', ha='center', va='bottom', fontsize=9)

            # График 2: P-value
            plt.subplot(1, 2, 2)
            p_values = df['p_value'].tolist()
            log_p = [-np.log10(p) if p > 0 else 3 for p in p_values]

            plt.bar(range(len(constants)), log_p, color='steelblue', alpha=0.7)
            plt.axhline(y=-np.log10(0.05), color='red', linestyle='--', label='p=0.05')
            plt.axhline(y=-np.log10(0.01), color='darkred', linestyle='--', label='p=0.01')
            plt.xticks(range(len(constants)), constants, rotation=45, ha='right')
            plt.ylabel('-log10(p-value)')
            plt.title('Статистическая значимость')
            plt.legend()

            plt.tight_layout()
            plt.savefig(output_file, dpi=150, bbox_inches='tight')
            plt.show()

            print(f"📊 График сохранён: {output_file}")

        except Exception as e:
            print(f"⚠ Ошибка при построении графика: {e}")

def main():
    """
    Основной пайплайн анализа.
    """
    print("🧬 GENOMIC CONSTANTS MAPPING TOOL v2.1")
    print("="*70)
    print("Автор: Андрей Шимельфениг")
    print("Дата: 2026")
    print("="*70)

    # Инициализация
    fetcher = GenomicSequenceFetcher(cache_dir="./genome_cache")
    matcher = ConstantsMatcher(fetcher)

    # Константы для анализа
    constants = [
        PhysicalConstant("Скорость света (c)", "299792458",
                        "ASPA", "chr10", 100188400, 100188418,
                        "AGGCGCCTGCAGCACCGA"),
        PhysicalConstant("Постоянная Планка (h)", "6.626070",
                        "MT-ND5", "chrM", 12337, 12351,
                        "CGCGAGCGAACTAA"),
        PhysicalConstant("Масса электрона (m_e)", "9.10938",
                        "SCN9A", "chr2", 166210400, 166210413,
                        "GCACAAAGCATGA"),
        PhysicalConstant("Гравитационная постоянная (G)", "6.67430",
                        "COL1A1", "chr17", 50183200, 50183212,
                        "CGCGCTCAATAA"),
        PhysicalConstant("Постоянная Больцмана (k)", "1.38064",
                        "TRPV1", "chr17", 3642010, 3642022,
                        "ACATGAACGACG"),
        PhysicalConstant("Постоянная тонкой структуры (α)", "137.035",
                        "KCNH2", "chr7", 150950100, 150950112,
                        "ACATCTAAATCC"),
        PhysicalConstant("Число Авогадро (N_A)", "6.02214",
                        "HBA1", "chr16", 176680, 176692,
                        "CGAAAGAGACCA")
    ]

    # Анализ каждой константы
    print("\n🔬 НАЧАЛО АНАЛИЗА")
    print("="*70)

    for const in constants:
        matcher.analyze_constant(const, debug=True)

    # Генерация отчёта
    print("\n" + "="*70)
    df = matcher.generate_report('genomic_constants_results')

    # Визуализация
    if len(df) > 0:
        matcher.plot_results(df, 'genomic_constants_plot.png')

    print("\n✨ Анализ завершён!")

    # Дополнительная информация
    print("\n📋 РЕКОМЕНДАЦИИ:")
    print("  1. Для повышения точности увеличьте число перестановок (сейчас 1000)")
    print("  2. Для полногеномного анализа используйте дополнительные гены")
    print("  3. Для публикации проведите анализ на 10,000+ перестановках")

    return df

if __name__ == "__main__":
    df_results = main()








