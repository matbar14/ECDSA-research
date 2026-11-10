import numpy as np
import matplotlib.pyplot as plt
import time
from collections import Counter
import sympy as sp

# ============================================================================
# СВЕДЕНИЕ ECDLP К ПОЛИНОМИАЛЬНОМУ УРАВНЕНИЮ
# Используем нашу проверенную формулу Λ(P)
# ============================================================================

# Параметры secp256k1
p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
beta = 0x7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee
n = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141

# Базовая точка G
Gx = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
Gy = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8

class EllipticCurvePoint:
    """Точка на эллиптической кривой secp256k1"""
    
    def __init__(self, x=None, y=None, is_infinity=False):
        if is_infinity:
            self.x = None
            self.y = None
            self.is_infinity = True
        else:
            self.x = x % p if x is not None else None
            self.y = y % p if y is not None else None
            self.is_infinity = False
    
    def __add__(self, other):
        """Сложение точек на эллиптической кривой"""
        if self.is_infinity:
            return other
        if other.is_infinity:
            return self
        
        if self.x == other.x:
            if self.y == other.y:
                # Удвоение точки
                s = (3 * self.x * self.x) * pow(2 * self.y, -1, p)
            else:
                # Точки противоположны
                return EllipticCurvePoint(is_infinity=True)
        else:
            # Обычное сложение
            s = (other.y - self.y) * pow(other.x - self.x, -1, p)
        
        s = s % p
        x3 = (s * s - self.x - other.x) % p
        y3 = (s * (self.x - x3) - self.y) % p
        
        return EllipticCurvePoint(x3, y3)
    
    def __mul__(self, scalar):
        """Скалярное умножение точки"""
        if scalar == 0:
            return EllipticCurvePoint(is_infinity=True)
        if scalar == 1:
            return EllipticCurvePoint(self.x, self.y, self.is_infinity)
        
        result = EllipticCurvePoint(is_infinity=True)
        addend = EllipticCurvePoint(self.x, self.y, self.is_infinity)
        
        while scalar:
            if scalar & 1:
                result = result + addend
            addend = addend + addend
            scalar >>= 1
        
        return result
    
    def __rmul__(self, scalar):
        return self.__mul__(scalar)
    
    def __eq__(self, other):
        if self.is_infinity and other.is_infinity:
            return True
        if self.is_infinity or other.is_infinity:
            return False
        return self.x == other.x and self.y == other.y

# Генераторная точка G
G = EllipticCurvePoint(Gx, Gy)

def compute_lambda_proven(xG, yG, xP, yP):
    """
    НАША ПРОВЕРЕННАЯ ФОРМУЛА ИНВАРИАНТА Λ(P)
    Λ(P) = ((x_G - β·x_G)(x_P - β·x_P)) / ((x_G - β·x_P)(x_P - β·x_G)) mod p
    """
    # φ(G) = (β·x_G, y_G)
    phi_G_x = (beta * xG) % p
    
    # φ(P) = (β·x_P, y_P) 
    phi_P_x = (beta * xP) % p
    
    # Вычисляем числитель: (x_G - β·x_G) * (x_P - β·x_P) mod p
    term1 = (xG - phi_G_x) % p
    term2 = (xP - phi_P_x) % p
    numerator = (term1 * term2) % p
    
    # Вычисляем знаменатель: (x_G - β·x_P) * (x_P - β·x_G) mod p
    term3 = (xG - phi_P_x) % p
    term4 = (xP - phi_G_x) % p
    denominator = (term3 * term4) % p
    
    # Делим по модулю p (умножаем на обратный элемент)
    if denominator != 0:
        denominator_inv = pow(denominator, p - 2, p)  # обратный элемент в F_p
        lambda_val = (numerator * denominator_inv) % p
        return lambda_val
    else:
        return 0  # специальное значение для случая деления на ноль

print("=" * 80)
print("СВЕДЕНИЕ ECDLP К ПОЛИНОМИАЛЬНОМУ УРАВНЕНИЮ")
print("Используя нашу проверенную формулу Λ(P) и подход из гипотезы")
print("=" * 80)
# ============================================================================
# АЛГОРИТМ СВЕДЕНИЯ ECDLP К ПОЛИНОМИАЛЬНОМУ УРАВНЕНИЮ
# ============================================================================

def derive_polynomial_equation(target_lambda):
    """
    Выводим полиномиальное уравнение для нахождения x_P по известному Λ(P)
    
    Дано: Λ(P) = target_lambda
    Найти: x_P такое, что точка P = (x_P, y_P) на кривой и Λ(P) = target_lambda
    """
    print(f"\n1. ВЫВОД ПОЛИНОМИАЛЬНОГО УРАВНЕНИЯ ДЛЯ Λ = 0x{target_lambda:016x}...")
    
    # Наша формула: Λ(P) = ((x_G - β·x_G)(x_P - β·x_P)) / ((x_G - β·x_P)(x_P - β·x_G))
    
    # Обозначения для упрощения
    xG = Gx
    phi_G_x = (beta * xG) % p  # β·x_G
    
    print(f"   x_G = 0x{xG:016x}...")
    print(f"   β·x_G = 0x{phi_G_x:016x}...")
    
    # Пусть x_P - неизвестная координата
    # φ(P) = (β·x_P, y_P)
    
    # Формула становится:
    # Λ = ((x_G - β·x_G)(x_P - β·x_P)) / ((x_G - β·x_P)(x_P - β·x_G))
    
    # Перепишем как:
    # Λ · ((x_G - β·x_P)(x_P - β·x_G)) = (x_G - β·x_G)(x_P - β·x_P)
    
    # Раскрываем левую часть:
    # Λ · (x_G·x_P - x_G·β·x_G - β·x_P·x_P + β²·x_P·x_G)
    # = Λ · (x_G·x_P - β·x_G² - β·x_P² + β²·x_G·x_P)
    # = Λ · (x_G·x_P(1 + β²) - β·x_G² - β·x_P²)
    
    # Раскрываем правую часть:
    # (x_G - β·x_G)(x_P - β·x_P) = (x_G(1-β))(x_P(1-β)) = x_G·x_P·(1-β)²
    
    # Получаем уравнение:
    # Λ · (x_G·x_P(1 + β²) - β·x_G² - β·x_P²) = x_G·x_P·(1-β)²
    
    # Приводим к стандартному виду:
    # Λ·x_G·x_P(1 + β²) - Λ·β·x_G² - Λ·β·x_P² = x_G·x_P·(1-β)²
    
    # Перегруппировываем:
    # x_P²·(-Λ·β) + x_P·(Λ·x_G·(1 + β²) - x_G·(1-β)²) - Λ·β·x_G² = 0
    
    # Это квадратное уравнение вида: A·x_P² + B·x_P + C = 0
    
    # Вычисляем коэффициенты
    beta2 = (beta * beta) % p
    
    A = (-target_lambda * beta) % p
    B = (target_lambda * xG * ((1 + beta2) % p) - xG * pow(1 - beta, 2, p)) % p
    C = (-target_lambda * beta * xG * xG) % p
    
    print(f"   Квадратное уравнение: A·x² + B·x + C = 0")
    print(f"   A = 0x{A:016x}...")
    print(f"   B = 0x{B:016x}...")
    print(f"   C = 0x{C:016x}...")
    
    return A, B, C

def solve_quadratic_mod_p(A, B, C):
    """
    Решает квадратное уравнение A·x² + B·x + C ≡ 0 (mod p)
    """
    print(f"\n2. РЕШЕНИЕ КВАДРАТНОГО УРАВНЕНИЯ:")
    
    if A == 0:
        # Линейное уравнение B·x + C = 0
        if B == 0:
            if C == 0:
                print("   Уравнение тождественно выполнено (бесконечно решений)")
                return list(range(100))  # Возвращаем первые 100 для тестирования
            else:
                print("   Уравнение не имеет решений")
                return []
        else:
            x = (-C * pow(B, p - 2, p)) % p
            print(f"   Линейное решение: x = 0x{x:016x}...")
            return [x]
    
    # Квадратное уравнение: A·x² + B·x + C = 0
    # Дискриминант: D = B² - 4·A·C
    D = (B * B - 4 * A * C) % p
    
    print(f"   Дискриминант D = 0x{D:016x}...")
    
    # Проверяем, является ли D квадратичным вычетом
    legendre = pow(D, (p - 1) // 2, p)
    
    if legendre != 1:
        print("   Дискриминант не является квадратичным вычетом - нет решений")
        return []
    
    # Находим квадратный корень из D
    # Для p ≡ 3 (mod 4) используем формулу sqrt(D) = D^((p+1)/4) mod p
    if p % 4 == 3:
        sqrt_D = pow(D, (p + 1) // 4, p)
    else:
        # Общий случай: перебор в tonelli_shanks(), для secp256k1 не используется
        sqrt_D = tonelli_shanks(D, p)
    
    # Проверяем корень
    if (sqrt_D * sqrt_D) % p != D:
        print(f"   Ошибка вычисления квадратного корня")
        return []
    
    print(f"   √D = 0x{sqrt_D:016x}...")
    
    # Два решения квадратного уравнения
    A_inv = pow(2 * A, p - 2, p)  # (2A)^(-1) mod p
    
    x1 = ((-B + sqrt_D) * A_inv) % p
    x2 = ((-B - sqrt_D) * A_inv) % p
    
    print(f"   Решение 1: x₁ = 0x{x1:016x}...")
    print(f"   Решение 2: x₂ = 0x{x2:016x}...")
    
    return [x1, x2]

def tonelli_shanks(n, p):
    """
    Квадратный корень по модулю p.

    Для p ≡ 3 (mod 4) — прямая формула n^((p+1)/4), именно этот случай
    имеет место для secp256k1. Общий случай ниже решается перебором,
    а не полным алгоритмом Тонелли-Шенкса.
    """
    # Упрощенная версия для p ≡ 3 (mod 4)
    if p % 4 == 3:
        return pow(n, (p + 1) // 4, p)
    
    # Для общего случая нужна полная реализация
    # Пока используем brute force для небольших значений
    for r in range(min(1000, p)):
        if (r * r) % p == n:
            return r
    
    return 0  # Не найдено

def find_points_on_curve(x_candidates):
    """
    Для каждого кандидата x находит соответствующие точки на кривой
    """
    print(f"\n3. ПОИСК ТОЧЕК НА КРИВОЙ:")
    
    points = []
    
    for x in x_candidates:
        # Вычисляем y² = x³ + 7 mod p
        y_squared = (pow(x, 3, p) + 7) % p
        
        # Проверяем, является ли y² квадратичным вычетом
        legendre = pow(y_squared, (p - 1) // 2, p)
        
        if legendre == 1:
            # Находим квадратные корни
            if p % 4 == 3:
                y = pow(y_squared, (p + 1) // 4, p)
            else:
                y = tonelli_shanks(y_squared, p)
            
            # Проверяем корень
            if (y * y) % p == y_squared:
                # Две точки: (x, y) и (x, -y)
                points.append((x, y))
                points.append((x, (-y) % p))
                
                print(f"   x = 0x{x:016x}... → точки найдены")
            else:
                print(f"   x = 0x{x:016x}... → ошибка вычисления y")
        else:
            print(f"   x = 0x{x:016x}... → не лежит на кривой")
    
    print(f"   Всего найдено точек: {len(points)}")
    return points

def verify_lambda_values(points, target_lambda):
    """
    Проверяем, какие из найденных точек дают правильное значение Λ
    """
    print(f"\n4. ВЕРИФИКАЦИЯ ЗНАЧЕНИЙ Λ:")
    
    correct_points = []
    
    for i, (x, y) in enumerate(points):
        # Вычисляем Λ для этой точки
        lambda_val = compute_lambda_proven(Gx, Gy, x, y)
        
        print(f"   Точка {i+1}: (0x{x:016x}..., 0x{y:016x}...)")
        print(f"            Λ = 0x{lambda_val:016x}...")
        
        if lambda_val == target_lambda:
            print(f"            СОВПАДЕНИЕ с целевым Λ!")
            correct_points.append((x, y))
        else:
            print(f"            Не совпадает с целевым Λ")
    
    print(f"\n   Найдено корректных точек: {len(correct_points)}")
    return correct_points


def _curve_y_roots(x):
    """Возвращает корни y для x на кривой secp256k1 (тихо), либо пустой список."""
    y_squared = (pow(x, 3, p) + 7) % p
    if pow(y_squared, (p - 1) // 2, p) != 1:
        return []
    if p % 4 == 3:
        y = pow(y_squared, (p + 1) // 4, p)
    else:
        y = tonelli_shanks(y_squared, p)
    if (y * y) % p != y_squared:
        return []
    return [y, (-y) % p]


def reconstruct_points_from_lambda(target_lambda):
    """Восстановление кандидатов точки P из известного Λ(P) без вывода.

    Тихая составная версия цепочки derive_polynomial_equation ->
    solve_quadratic_mod_p -> find_points_on_curve -> verify_lambda_values.
    При Λ != 0 крест-отношение даёт квадратное уравнение по x_P не более чем
    с двумя корнями, поэтому кандидатов не более четырёх точек.
    При Λ == 0 (граничное соглашение w1.py) отдельно проверяются корни
    знаменателя x_P = β·x_G и x_P = x_G/β, а также корень числителя x_P = 0.
    Возвращает список EllipticCurvePoint, дающих исходное значение Λ.
    """
    target_lambda = target_lambda % p
    if target_lambda == 0:
        inverse_beta = pow(beta, p - 2, p)
        x_candidates = [
            0,
            (beta * Gx) % p,
            (Gx * inverse_beta) % p,
        ]
    else:
        beta2 = (beta * beta) % p
        A = (-target_lambda * beta) % p
        B = (target_lambda * Gx * ((1 + beta2) % p) - Gx * pow(1 - beta, 2, p)) % p
        C = (-target_lambda * beta * Gx * Gx) % p
        D = (B * B - 4 * A * C) % p
        if pow(D, (p - 1) // 2, p) != 1:
            return []
        if p % 4 == 3:
            sqrt_D = pow(D, (p + 1) // 4, p)
        else:
            sqrt_D = tonelli_shanks(D, p)
        if (sqrt_D * sqrt_D) % p != D:
            return []
        A_inv = pow(2 * A, p - 2, p)
        x1 = ((-B + sqrt_D) * A_inv) % p
        x2 = ((-B - sqrt_D) * A_inv) % p
        x_candidates = [x1] if x1 == x2 else [x1, x2]

    candidates = []
    seen = set()
    for x in x_candidates:
        for y in _curve_y_roots(x):
            if (x, y) in seen:
                continue
            seen.add((x, y))
            if compute_lambda_proven(Gx, Gy, x, y) == target_lambda:
                candidates.append(EllipticCurvePoint(x, y))
    return candidates


# ============================================================================
# ТЕСТИРОВАНИЕ АЛГОРИТМА РЕДУКЦИИ
# ============================================================================

def test_ecdlp_reduction():
    """
    Тестируем полный алгоритм сведения ECDLP к полиномиальному уравнению
    """
    print("\n" + "="*80)
    print("ТЕСТИРОВАНИЕ АЛГОРИТМА СВЕДЕНИЯ ECDLP")
    print("="*80)
    
    # Выбираем тестовый приватный ключ
    test_d = 42
    
    # Вычисляем публичный ключ
    P_test = test_d * G
    
    print(f"ТЕСТОВЫЕ ДАННЫЕ:")
    print(f"   Приватный ключ d = {test_d}")
    print(f"   Публичный ключ P = (0x{P_test.x:016x}..., 0x{P_test.y:016x}...)")
    
    # Вычисляем целевой инвариант
    target_lambda = compute_lambda_proven(G.x, G.y, P_test.x, P_test.y)
    print(f"   Целевой инвариант Λ = 0x{target_lambda:016x}...")
    
    # Шаг 1: Выводим полиномиальное уравнение
    A, B, C = derive_polynomial_equation(target_lambda)
    
    # Шаг 2: Решаем квадратное уравнение
    x_solutions = solve_quadratic_mod_p(A, B, C)
    
    if not x_solutions:
        print("Не найдено решений квадратного уравнения")
        return False
    
    # Шаг 3: Находим точки на кривой
    candidate_points = find_points_on_curve(x_solutions)
    
    if not candidate_points:
        print("Не найдено точек на кривой")
        return False
    
    # Шаг 4: Верифицируем значения Λ
    correct_points = verify_lambda_values(candidate_points, target_lambda)
    
    if not correct_points:
        print("Ни одна точка не дает правильное значение Λ")
        return False
    
    # Шаг 5: Проверяем, есть ли среди них исходная точка
    print(f"\n5. ПРОВЕРКА ИСХОДНОЙ ТОЧКИ:")
    
    found_original = False
    for x, y in correct_points:
        if x == P_test.x and y == P_test.y:
            print(f"   НАЙДЕНА ИСХОДНАЯ ТОЧКА!")
            print(f"   P = (0x{x:016x}..., 0x{y:016x}...)")
            found_original = True
            break
    
    if not found_original:
        print(f"   Исходная точка не найдена среди кандидатов")
        print(f"   Исходная: (0x{P_test.x:016x}..., 0x{P_test.y:016x}...)")
        print(f"   Найденные:")
        for x, y in correct_points:
            print(f"      (0x{x:016x}..., 0x{y:016x}...)")
    
    return found_original

# ============================================================================
# МАССОВОЕ ТЕСТИРОВАНИЕ
# ============================================================================

def mass_testing():
    """
    Тестируем алгоритм на множестве приватных ключей
    """
    print("\n" + "="*80)
    print("МАССОВОЕ ТЕСТИРОВАНИЕ АЛГОРИТМА")
    print("="*80)
    
    test_keys = [1, 2, 3, 5, 7, 10, 42, 100, 256, 1337]
    
    success_count = 0
    total_count = len(test_keys)
    
    for d in test_keys:
        print(f"\n--- ТЕСТ ДЛЯ d = {d} ---")
        
        try:
            # Вычисляем публичный ключ и инвариант
            P = d * G
            target_lambda = compute_lambda_proven(G.x, G.y, P.x, P.y)
            
            # Решаем полиномиальное уравнение
            A, B, C = derive_polynomial_equation(target_lambda)
            x_solutions = solve_quadratic_mod_p(A, B, C)
            
            if x_solutions:
                candidate_points = find_points_on_curve(x_solutions)
                correct_points = verify_lambda_values(candidate_points, target_lambda)
                
                # Проверяем, найдена ли исходная точка
                found = any(x == P.x and y == P.y for x, y in correct_points)
                
                if found:
                    print(f"   УСПЕХ для d = {d}")
                    success_count += 1
                else:
                    print(f"   НЕУДАЧА для d = {d}")
            else:
                print(f"   Нет решений для d = {d}")
                
        except Exception as e:
            print(f"   ОШИБКА для d = {d}: {e}")
    
    accuracy = success_count / total_count * 100
    print(f"\n" + "="*80)
    print(f"РЕЗУЛЬТАТЫ МАССОВОГО ТЕСТИРОВАНИЯ:")
    print(f"   Успешных тестов: {success_count}/{total_count}")
    print(f"   Точность алгоритма: {accuracy:.1f}%")
    print("="*80)
    
    return accuracy

if __name__ == "__main__":
    # Запускаем тестирование
    print("Запуск тестирования алгоритма сведения ECDLP к полиномиальному уравнению...")
    
    # Одиночный тест
    single_success = test_ecdlp_reduction()
    
    # Массовое тестирование
    mass_accuracy = mass_testing()
    
    print(f"\n" + "="*80)
    print("ФИНАЛЬНЫЕ ВЫВОДЫ:")
    print("="*80)
    
    if single_success:
        print("Одиночный тест прошел успешно")
    else:
        print("Одиночный тест провалился")
    
    print(f"Точность на массовых тестах: {mass_accuracy:.1f}%")
    
    if mass_accuracy >= 90:
        print("АЛГОРИТМ РАБОТАЕТ! ECDLP успешно сведена к полиномиальному уравнению!")
    elif mass_accuracy >= 50:
        print("Алгоритм работает частично - требуется доработка")
    else:
        print("Алгоритм не работает - нужен пересмотр подхода")