from itertools import product
from collections import deque


class PerformListAlignment:

    def __init__(self, x, y):
        """
        Initialize
        :param x: list 1
        :param y: list 2
        :return:
        """
        self.x = x
        self.y = y
        self.upperList = []
        self.lowerList = []
        self.matchScore = 0
        self.mismatchScore = 3
        self.gapScore = 7
        self.totalScore = 0

    def needleman_wunsch(self):
        """
        Run the Needleman-Wunsch algorithm on two sequences.
        Code based on pseudocode in Section 3 of:
        :param x:
        :param y:
        :return:
        """

        N, M = len(self.x), len(self.y)
        s = lambda a, b: int(a == b)

        DIAG = -1, -1
        LEFT = -1, 0
        UP = 0, -1

        # Create tables F and Ptr
        F = {}
        Ptr = {}

        F[-1, -1] = 0
        for i in range(N):
            F[i, -1] = -i
        for j in range(M):
            F[-1, j] = -j

        option_Ptr = DIAG, LEFT, UP
        for i, j in product(range(N), range(M)):
            option_F = (
                F[i - 1, j - 1] + s(self.x[i], self.y[j]),
                F[i - 1, j] - 1,
                F[i, j - 1] - 1,
            )
            F[i, j], Ptr[i, j] = max(zip(option_F, option_Ptr))

        # Work backwards from (N - 1, M - 1) to (0, 0)
        # to find the best alignment.
        alignment = deque()
        i, j = N - 1, M - 1
        while i >= 0 and j >= 0:
            direction = Ptr[i, j]
            if direction == DIAG:
                element = i, j
            elif direction == LEFT:
                element = i, None
            elif direction == UP:
                element = None, j
            alignment.appendleft(element)
            di, dj = direction
            i, j = i + di, j + dj
        while i >= 0:
            alignment.appendleft((i, None))
            i -= 1
        while j >= 0:
            alignment.appendleft((None, j))
            j -= 1
        return list(alignment)

    # def print_alignment(self, x, y, alignment):
    #     print("*".join(
    #         "-" if i is None else x[i] for i, _ in alignment
    #     ))
    #     print("*".join(
    #         "-" if j is None else y[j] for _, j in alignment
    #     ))
    def calculateScore(self):
        score = 0
        for i,u in enumerate(self.upperList):
            if self.upperList[i] is None or self.lowerList[i] is None:
                # gap
                score += self.gapScore
            elif self.upperList[i] == self.lowerList[i]:
                # matched
                score += self.matchScore
            else:
                # mismatched
                score += self.mismatchScore
        self.totalScore = score

    def run(self):
        alignment = self.needleman_wunsch()
        self.upperList = [None if i is None else self.x[i] for i, _ in alignment]
        self.lowerList = [None if j is None else self.y[j] for _, j in alignment]
        self.calculateScore()
        # print("*".join("-" if i is None else self.x[i] for i, _ in alignment))
        # print("*".join("-" if j is None else self.y[j] for _, j in alignment))

    # def align_fast(self, x, y):
    #     """Align two sequences, maximizing the
    #     alignment score, using the Needleman-Wunsch
    #     algorithm.
    #
    #     x, y -- sequences.
    #     """
    #     return needleman_wunsch(x, y)


# a = "CAT"
# b = "CT"

# a = ["C", "A", "T"]
# b = ["C", "T"]

# a = ["dog", "cat", "bat"]
# b = ["dog", "bat"]

if __name__ == '__main__':
    a = ['dog', 'cat', 'monkey', 'mouse', 'elephant','axa','eiei']
    b = ['cat', 'penguin', 'elephant']
    obj = PerformListAlignment(a, b)
    obj.run()
    # print_alignment(a, b, align_fast(a, b))
    print(obj.upperList)
    print(obj.lowerList)
    print(obj.totalScore)
