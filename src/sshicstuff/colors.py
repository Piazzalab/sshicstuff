import random



def generate(type: str, n: int, a: float = 0.8, seed: int = 42) -> list[str]:
    random.seed(seed)

    if type == 'hex':
        white = '#FFFFFF'
        blue = '#0080FF'
        res = [f"#{random.randint(0, 0xFFFFFF):06x}" for _ in range(n)]
        if white in res:
            res[res.index(white)] = blue

        return res

    elif type == 'rgba':
        white = 'rgba(255, 255, 255, 0.8)'
        blue = 'rgba(0, 128, 255, 0.8)'

        res = [
            (f"rgba({random.randint(0, 255)}, "
             f"{random.randint(0, 255)}, "
             f"{random.randint(0, 255)}, {a})") for _ in range(n)
        ]

        if white in res:
            res[res.index(white)] = blue

        return res

    else:
        raise ValueError("type must be 'hex' or 'rgba'")




chr_colorbar = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
    '#bcbd22', '#17becf', '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5', '#c49c94',
    '#f7b6d2', '#c7c7c7', '#dbdb8d', '#9edae5', '#393b79', '#5254a3', '#6b6ecf', '#9c9ede',
]