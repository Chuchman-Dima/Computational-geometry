import tkinter as tk
from tkinter import ttk, messagebox
import random
import time

# --- ВАЖЛИВО: Потрібні бібліотеки! ---
try:
    import numpy as np
    from scipy.spatial import Voronoi, Delaunay, QhullError
except ImportError:
    messagebox.showerror("Помилка імпорту",
                         "Необхідні бібліотеки 'scipy' та 'numpy' не знайдено.\n"
                         "Будь ласка, встановіть їх, виконавши команду:\n"
                         "pip install scipy numpy")
    exit()


def findLER_Algorithm_1_Naive(points, bounds):
    """
    АЛГОРИТМ 1: Демонстрація "Методу Орловського"
    """
    x_coords = sorted(list(set([bounds[0]] + [p[0] for p in points] + [bounds[1]])))
    y_coords = sorted(list(set([bounds[2]] + [p[1] for p in points] + [bounds[3]])))

    n_x = len(x_coords)
    n_y = len(y_coords)
    maxArea = 0
    bestRect = []

    for i in range(n_x - 1):
        for j in range(i + 1, n_x):
            for k in range(n_y - 1):
                for l in range(k + 1, n_y):
                    x_min, x_max = x_coords[i], x_coords[j]
                    y_min, y_max = y_coords[k], y_coords[l]

                    is_empty = True
                    for p in points:
                        if (p[0] > x_min and p[0] < x_max and \
                                p[1] > y_min and p[1] < y_max):
                            is_empty = False
                            break

                    if is_empty:
                        area = (x_max - x_min) * (y_max - y_min)
                        if area > maxArea:
                            maxArea = area
                            bestRect = [x_min, y_min, x_max, y_max]

    return bestRect, maxArea, None


def findLER_Algorithm_2_Voronoi_Demo(points, bounds):
    """
    АЛГОРИТМ 2: Демонстрація "Методу Діаграми Вороного"
    """

    # --- Обчислення діаграми Вороного для візуалізації ---
    details = {'type': 'voronoi', 'data': None}
    try:
        if len(points) >= 3:
            np_points = np.array(points)
            width = bounds[1] - bounds[0]
            height = bounds[3] - bounds[2]
            scale = max(width, height) * 50

            ghost_points = np.array([
                [bounds[0] - scale, bounds[2] - scale],
                [bounds[1] + scale, bounds[2] - scale],
                [bounds[0] - scale, bounds[3] + scale],
                [bounds[1] + scale, bounds[3] + scale]
            ])

            full_point_set = np.concatenate((np_points, ghost_points))
            vor = Voronoi(full_point_set)
            details['data'] = vor
        else:
            details = None
    except (QhullError, ValueError) as e:
        print(f"Помилка при побудові діаграми Вороного: {e}")
        details = None

    y_candidates = sorted(list(set([bounds[2]] + [p[1] for p in points] + [bounds[3]])))
    n_y = len(y_candidates)

    maxArea = 0
    bestRect = []

    for i in range(n_y - 1):
        for j in range(i + 1, n_y):
            y_min, y_max = y_candidates[i], y_candidates[j]
            height = y_max - y_min

            strip_points_x = [p[0] for p in points if p[1] > y_min and p[1] < y_max]
            x_in_strip = sorted(list(set([bounds[0]] + strip_points_x + [bounds[1]])))

            max_width = 0
            x_left = bounds[0]

            for k in range(len(x_in_strip) - 1):
                width = x_in_strip[k + 1] - x_in_strip[k]
                if width > max_width:
                    max_width = width
                    x_left = x_in_strip[k]

            area = max_width * height

            if area > maxArea:
                maxArea = area
                bestRect = [x_left, y_min, x_left + max_width, y_max]

    return bestRect, maxArea, details


def findLER_Algorithm_3_Triangulation_Demo(points, bounds):
    """
    АЛГОРИТМ 3: Демонстрація "Методу Тріангуляції"
    """

    # --- Обчислення тріангуляції для візуалізації ---
    details = {'type': 'triangulation', 'data': None}
    if len(points) >= 3:
        try:
            tri = Delaunay(np.array(points))
            details['data'] = tri
        except (QhullError, ValueError) as e:
            print(f"Помилка при побудові тріангуляції: {e}")
            details = None

    x_coords = sorted(list(set([bounds[0]] + [p[0] for p in points] + [bounds[1]])))
    n_x = len(x_coords)
    maxArea = 0
    bestRect = []

    for i in range(n_x - 1):
        for j in range(i + 1, n_x):
            x_min = x_coords[i]
            x_max = x_coords[j]
            width = x_max - x_min

            strip_points_y = [p[1] for p in points if p[0] > x_min and p[0] < x_max]

            y_in_strip = sorted(list(set([bounds[2]] + strip_points_y + [bounds[3]])))

            max_height = 0
            y_bottom = bounds[2]

            for k in range(len(y_in_strip) - 1):
                height = y_in_strip[k + 1] - y_in_strip[k]
                if height > max_height:
                    max_height = height
                    y_bottom = y_in_strip[k]

            area = width * max_height
            if area > maxArea:
                maxArea = area
                bestRect = [x_min, y_bottom, x_max, y_bottom + max_height]

    return bestRect, maxArea, details


# --- РОЗДІЛ: ГРАФІЧНИЙ ІНТЕРФЕЙС (GUI) ---

class LerApplication(tk.Tk):

    def __init__(self):
        super().__init__()
        self.title("Демонстрація: Найбільший Порожній Прямокутник (РМП (АО6))")
        self.geometry("1000x700")

        self.points_list = []
        self.bounds = [10, 690, 10, 690]
        self.np_points = np.array([])

        self.create_widgets()

    def create_widgets(self):
        control_frame = ttk.Frame(self, padding=10)
        control_frame.pack(side=tk.LEFT, fill=tk.Y)

        ttk.Label(control_frame, text="Керування", font=("Arial", 16, "bold")).pack(pady=10)

        input_frame = ttk.LabelFrame(control_frame, text="1. Режим вводу")
        input_frame.pack(fill=tk.X, padx=5, pady=5)
        self.input_mode = tk.StringVar(value="manual")
        ttk.Radiobutton(input_frame, text="Ручний (клік на полі) (до 100)", variable=self.input_mode,
                        value="manual").pack(anchor=tk.W)
        ttk.Radiobutton(input_frame, text="Автоматичний (до 10000)", variable=self.input_mode, value="auto").pack(
            anchor=tk.W)

        gen_frame = ttk.LabelFrame(control_frame, text="2. Генерація (для 'Авто')")
        gen_frame.pack(fill=tk.X, padx=5, pady=5)
        ttk.Label(gen_frame, text="Кількість точок:").pack(side=tk.LEFT, padx=5)
        self.num_points_var = tk.StringVar(value="50")
        ttk.Entry(gen_frame, textvariable=self.num_points_var, width=8).pack(side=tk.LEFT)
        ttk.Button(gen_frame, text="Згенерувати", command=self.generate_points).pack(pady=5, fill=tk.X)

        algo_frame = ttk.LabelFrame(control_frame, text="3. Вибір алгоритму")
        algo_frame.pack(fill=tk.X, padx=5, pady=5)
        self.algo_mode = tk.StringVar(value="voronoi")

        ttk.Radiobutton(algo_frame, text="Метод Орловського", variable=self.algo_mode,
                        value="orlowski").pack(anchor=tk.W)
        ttk.Radiobutton(algo_frame, text="Метод Тріангуляції", variable=self.algo_mode,
                        value="triangulation").pack(anchor=tk.W)
        ttk.Radiobutton(algo_frame, text="Метод Діаграми Вороного", variable=self.algo_mode,
                        value="voronoi").pack(anchor=tk.W)

        self.show_details_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            control_frame,
            text="Показувати деталі візуалізації\n(Тріангуляція / Діаграма Вороного)",
            variable=self.show_details_var
        ).pack(pady=5)

        ttk.Button(control_frame, text="РОЗРАХУВАТИ", style="Accent.TButton", command=self.run_algorithm).pack(pady=10,
                                                                                                               fill=tk.X)
        ttk.Button(control_frame, text="ОЧИСТИТИ", command=self.clear_canvas).pack(pady=5, fill=tk.X)

        self.style = ttk.Style()
        self.style.configure("Accent.TButton", foreground="black", background="yellow")

        self.status_var = tk.StringVar(value="Готово. Оберіть режим вводу.")
        status_bar = ttk.Label(self, textvariable=self.status_var, relief=tk.SUNKEN, anchor=tk.W, padding=5)
        status_bar.pack(side=tk.BOTTOM, fill=tk.X)

        self.canvas = tk.Canvas(self, bg="white", highlightthickness=1, highlightbackground="black")
        self.canvas.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=10, pady=10)

        self.canvas.bind("<Button-1>", self.on_canvas_click)
        self.canvas.bind("<Configure>", self.update_bounds)
        self.after(100, self.initial_draw)

    def initial_draw(self):
        self.update_bounds(
            type('Event', (object,), {'width': self.canvas.winfo_width(), 'height': self.canvas.winfo_height()}))

    def update_bounds(self, event):
        self.bounds = [10, event.width - 10, 10, event.height - 10]
        self.canvas.delete("bounds_rect")
        self.canvas.create_rectangle(
            self.bounds[0], self.bounds[2], self.bounds[1], self.bounds[3],
            outline="gray", dash=(2, 2), tags="bounds_rect"
        )
        self.canvas.tag_lower("bounds_rect")

    def draw_point(self, x, y):
        self.canvas.create_oval(x - 2.5, y - 2.5, x + 2.5, y + 2.5, fill="blue", outline="blue", tags="point")
        self.points_list.append((x, y))

    def on_canvas_click(self, event):
        if self.input_mode.get() == "manual":
            if len(self.points_list) >= 100:
                self.status_var.set("Помилка: В ручному режимі до 100 точок.")
                return
            x, y = event.x, event.y
            if (self.bounds[0] < x < self.bounds[1] and
                    self.bounds[2] < y < self.bounds[3]):
                self.draw_point(x, y)
                self.np_points = np.array(self.points_list)
                self.status_var.set(f"Додано точку {len(self.points_list)}: ({x}, {y})")
            else:
                self.status_var.set("Клікніть всередині пунктирної рамки.")
        else:
            self.status_var.set("Перемкніться в 'Ручний режим' для додавання точок.")

    def generate_points(self):
        if self.input_mode.get() != "auto":
            self.status_var.set("Помилка: Спочатку оберіть 'Автоматичний' режим.")
            return
        try:
            n = int(self.num_points_var.get())
            if not (1 <= n <= 10000):
                raise ValueError
        except ValueError:
            messagebox.showerror("Помилка", "Кількість точок має бути числом від 1 до 10000.")
            return

        self.clear_canvas()
        self.status_var.set(f"Генерація {n} точок...")
        self.update_idletasks()

        b = self.bounds
        for _ in range(n):
            x = random.uniform(b[0], b[1])
            y = random.uniform(b[2], b[3])
            self.draw_point(x, y)

        self.np_points = np.array(self.points_list)
        self.status_var.set(f"Згенеровано {n} точок. Натисніть 'РОЗРАХУВАТИ'.")

    def clear_canvas(self):
        self.points_list = []
        self.np_points = np.array([])
        self.canvas.delete("all")
        self.update_bounds(
            type('Event', (object,), {'width': self.canvas.winfo_width(), 'height': self.canvas.winfo_height()}))
        self.status_var.set("Очищено. Можна починати.")

    def run_algorithm(self):
        if not self.points_list:
            self.status_var.set("Помилка: Немає точок для аналізу.")
            return

        if len(self.points_list) < 3 and self.algo_mode.get() != "orlowski":
            self.status_var.set("Помилка: Потрібно мінімум 3 точки для візуалізації.")
            return

        n = len(self.points_list)
        algo = self.algo_mode.get()

        if n > 100 and algo == "orlowski":
            if not messagebox.askyesno("Попередження",
                                       f"Ви запускаєте 'Орловського' з {n} точками. Це може зайняти ДУЖЕ багато часу. Продовжити?"): return
        if n > 500 and algo == "triangulation":
            if not messagebox.askyesno("Попередження",
                                       f"Ви запускаєте 'Тріангуляції' з {n} точками. Це може зайняти деякий час. Продовжити?"): return
        if n > 3000 and algo == "voronoi":
            if not messagebox.askyesno("Попередження",
                                       f"Ви запускаєте 'Вороного' з {n} точками. Це може зайняти деякий час. Продовжити?"): return

        self.status_var.set(f"Розрахунок... (n={n}, алгоритм: {algo})")
        self.update_idletasks()

        start_time = time.time()

        details = None
        if algo == "orlowski":
            bestRect, maxArea, details = findLER_Algorithm_1_Naive(self.points_list, self.bounds)
        elif algo == "triangulation":
            bestRect, maxArea, details = findLER_Algorithm_3_Triangulation_Demo(self.points_list, self.bounds)
        else:  # "voronoi_demo"
            # --- ВИПРАВЛЕНО ВИКЛИК ФУНКЦІЇ ---
            bestRect, maxArea, details = findLER_Algorithm_2_Voronoi_Demo(self.points_list, self.bounds)

        end_time = time.time()

        self.canvas.delete("result_rect")
        self.canvas.delete("details")

        if bestRect:
            x1, y1, x2, y2 = bestRect
            self.canvas.create_rectangle(
                x1, y1, x2, y2,
                outline="red", width=2, tags="result_rect"
            )

            if self.show_details_var.get() and details and details['data'] is not None:
                self.draw_details(details)

            self.status_var.set(f"Готово! Час: {end_time - start_time:.4f} c. Площа: {maxArea:.2f}")
        else:
            self.status_var.set(f"Готово! Час: {end_time - start_time:.4f} c. Не знайдено порожніх прямокутників.")

    def draw_details(self, details):
        """Малює візуалізацію (Діаграму або Тріангугуляцію)"""

        if details['type'] == 'triangulation':
            tri = details['data']
            for simplex in tri.simplices:
                p0, p1, p2 = self.np_points[simplex]
                self.canvas.create_line(p0[0], p0[1], p1[0], p1[1], fill="lightblue", tags="details")
                self.canvas.create_line(p1[0], p1[1], p2[0], p2[1], fill="lightblue", tags="details")
                self.canvas.create_line(p2[0], p2[1], p0[0], p0[1], fill="lightblue", tags="details")

        elif details['type'] == 'voronoi':
            vor = details['data']
            num_original_points = len(self.points_list)

            for (p1_idx, p2_idx), (v1_idx, v2_idx) in zip(vor.ridge_points, vor.ridge_vertices):

                # Не малюємо ребра, які з'єднують ДВІ "примарні" точки
                if p1_idx >= num_original_points and p2_idx >= num_original_points:
                    continue

                if v1_idx != -1 and v2_idx != -1:
                    v0 = vor.vertices[v1_idx]
                    v1 = vor.vertices[v2_idx]
                    self.canvas.create_line(
                        v0[0], v0[1], v1[0], v1[1],
                        fill="lightgreen", tags="details"
                    )

        self.canvas.tag_lower("details")  # Ховаємо під точки


# --- Запуск програми ---
if __name__ == "__main__":
    app = LerApplication()
    app.mainloop()