import tkinter as tk
import random

APP_TITLE = "Will you be my girlfriend?"

root = tk.Tk()
root.title("A Question")
root.geometry("700x420")
root.resizable(False, False)
root.configure(bg="#fff7f0")

header = tk.Label(
    root,
    text=APP_TITLE,
    font=("Segoe UI", 24, "bold"),
    bg="#fff7f0",
    fg="#c43d3d",
)
header.pack(pady=20)

subtext = tk.Label(
    root,
    text="Please say yes. No pressure though :) ",
    font=("Segoe UI", 12),
    bg="#fff7f0",
    fg="#5b4b4b",
)
subtext.pack()

button_area = tk.Frame(root, bg="#fff7f0", width=700, height=260)
button_area.pack(pady=20, fill="both", expand=True)

button_area.update_idletasks()

base_yes_size = 16
current_yes_size = base_yes_size
base_window_width = 700
base_window_height = 420

messages = [
    "Are you sure?",
    "Maybe try yes?",
    "Nope, still asking...",
    "Ok but the answer is yes.",
    "Come on, just click it.",
]


def move_no_button():
    button_area.update_idletasks()
    max_x = max(button_area.winfo_width() - 120, 20)
    max_y = max(button_area.winfo_height() - 60, 20)
    new_x = random.randint(20, max_x)
    new_y = random.randint(20, max_y)
    no_button.place(x=new_x, y=new_y)


def on_no():
    global current_yes_size
    current_yes_size += 4
    yes_button.config(font=("Segoe UI", current_yes_size, "bold"), padx=10, pady=8)
    ensure_window_fits_yes()
    see = messages[current_yes_size % len(messages)]
    subtext.config(text=see)
    move_no_button()


def on_yes():
    header.config(text="Yay! I knew it! ")
    subtext.config(text="You made me the happiest person ever.")
    yes_button.config(state="disabled")
    no_button.config(state="disabled")

def ensure_window_fits_yes():
    root.update_idletasks()
    yes_w = yes_button.winfo_reqwidth()
    yes_h = yes_button.winfo_reqheight()
    area_w = button_area.winfo_width()
    area_h = button_area.winfo_height()

    extra_w = max(0, yes_w + 80 - area_w)
    extra_h = max(0, yes_h + 80 - area_h)

    if extra_w or extra_h:
        new_w = max(root.winfo_width() + extra_w, base_window_width)
        new_h = max(root.winfo_height() + extra_h, base_window_height)
        root.geometry(f"{new_w}x{new_h}")


yes_button = tk.Button(
    button_area,
    text="Yes",
    font=("Segoe UI", base_yes_size, "bold"),
    bg="#ff6b6b",
    fg="white",
    activebackground="#ff8787",
    activeforeground="white",
    relief="flat",
    command=on_yes,
    padx=10,
    pady=6,
)

yes_button.place(x=210, y=110)

no_button = tk.Button(
    button_area,
    text="No",
    font=("Segoe UI", 14),
    bg="#f4f4f4",
    fg="#333333",
    activebackground="#ededed",
    relief="flat",
    command=on_no,
    padx=8,
    pady=4,
)

no_button.place(x=390, y=120)

root.mainloop()
