import os
from dash_iconify import DashIconify
import dash_mantine_components as dmc

#   source :
#   https://community.plotly.com/t/file-explorer-tree-generator-for-local-files/68732


class FileTree:

    def __init__(self, filepath: os.PathLike):
        """
        Usage: component = FileTree('Path/to/my/File').render()
        """
        self.filepath = filepath

    def render(self) -> dmc.Accordion:
        return dmc.Accordion(
            self.build_tree(self.filepath, isRoot=True),
            multiple=True)

    @staticmethod
    def flatten(self, l):
        return [item for sublist in l for item in sublist]

    @staticmethod
    def make_file(self, file_name):
        return dmc.Text(
            [DashIconify(icon="akar-icons:file"), " ", file_name], style={"paddingTop": '5px'}
        )

    @staticmethod
    def make_folder(self, folder_name):
        return [DashIconify(icon="akar-icons:folder"), " ", folder_name]

    def build_tree(self, path, isRoot=False):
        d = []
        if os.path.isdir(path):
            children = self.flatten([self.build_tree(os.path.join(path, x)) for x in os.listdir(path)])
            if isRoot:
                d.append(
                    dmc.AccordionItem(
                        children=children,
                        label=self.make_folder(os.path.basename(path)))
                )
            else:
                d.append(
                    dmc.Accordion(children=[
                        dmc.AccordionItem(
                            children=children,
                            label=self.make_folder(os.path.basename(path)))
                    ],
                        multiple=True)
                )
        else:
            d.append(self.make_file(os.path.basename(path)))
        return d
