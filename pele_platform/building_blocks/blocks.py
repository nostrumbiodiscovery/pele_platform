from abc import ABC, abstractmethod


class Block(ABC):
    @abstractmethod
    def run(self):
        pass
