#include <gtk/gtk.h>

class VisSimParams;

void on_button1_clicked(GtkButton * widget, VisSimParams * p);

void on_button2_clicked(GtkButton * widget);

void on_targetchooser_file_set(GtkFileChooser * widget, VisSimParams * p);

void on_arraychooser_file_set(GtkFileChooser * widget, VisSimParams * p);

void on_instrumentchooser_file_set(GtkFileChooser * widget, VisSimParams * p);

void on_SMchooser_file_set(GtkFileChooser * widget, VisSimParams * p);

void on_HAchooser_file_set(GtkFileChooser * widget, VisSimParams * p);

void on_oifitschooser_changed(GtkEntry * widget, VisSimParams * p);

void gui_main(int argc, char *argv[]);
