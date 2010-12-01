#include <cstdlib>
#include <string>
#include <iostream>
#include <stdexcept>
#include <gtk/gtk.h>

#include "userinterface.h"
#include "VisSimParams.h"
#include "oifits_sim.h"

using std::string;
using std::cout;
using std::endl;

void on_button1_clicked(GtkButton * widget, VisSimParams * p)
{
	if (p->have_all_params())
	{
		cout << "Generating OIFITS from:" << endl;
		cout << "Target selected: " << p->target_filename << endl;
		cout << "Array selected: " << p->array_filename << endl;
		cout << "Spectral Mode selected: " << p->
		   SpectralMode_filename << endl;
		cout << "Instrument selected: " << p->instrument_filename << endl;
		cout << "Hour angles selected:" << p->HourAngles_filename << endl;
		cout << "The output file is " << p->oifits_filename << endl;
		try
		{
			run_sim(p);
		}
		catch(const std::runtime_error & ex)
		{
			GtkWidget *dialog = gtk_message_dialog_new(NULL,
			   GTK_DIALOG_DESTROY_WITH_PARENT,
			   GTK_MESSAGE_ERROR,
			   GTK_BUTTONS_CLOSE,
			   ex.what(),
			   NULL);
			gtk_dialog_run(GTK_DIALOG(dialog));
			gtk_widget_destroy(dialog);
		}
	}
	else
	{
		cout << "Filename(s) not specified" << endl;
	}
}

void on_button2_clicked(GtkButton * widget)
{
	gtk_main_quit();
}

void on_targetchooser_file_set(GtkFileChooser * widget, VisSimParams * p)
{
	p->target_filename = gtk_file_chooser_get_filename(widget);
}

void on_arraychooser_file_set(GtkFileChooser * widget, VisSimParams * p)
{
	p->array_filename = gtk_file_chooser_get_filename(widget);
}

void on_instrumentchooser_file_set(GtkFileChooser * widget, VisSimParams * p)
{
	p->instrument_filename = gtk_file_chooser_get_filename(widget);

}

void on_SMchooser_file_set(GtkFileChooser * widget, VisSimParams * p)
{
	p->SpectralMode_filename = gtk_file_chooser_get_filename(widget);

}

void on_HAchooser_file_set(GtkFileChooser * widget, VisSimParams * p)
{
	p->HourAngles_filename = gtk_file_chooser_get_filename(widget);

}

void on_oifitschooser_changed(GtkEntry * widget, VisSimParams * p)
{
	p->oifits_filename = gtk_entry_get_text(widget);
}

void gui_main(int argc, char *argv[])
{
	GtkBuilder *main_window = NULL;
	GError *error = NULL;
	GtkWidget *widget;
	VisSimParams p = VisSimParams();
	const char *gladeName = "GUI.glade";

	gtk_init(&argc, &argv);

	/*
	 * load the interface 
	 */
	if (g_file_test(gladeName, G_FILE_TEST_EXISTS))
	{
		main_window = gtk_builder_new();
		if (!gtk_builder_add_from_file(main_window, gladeName, &error))
		{
			g_warning("Couldn't load builder file: %s", error->message);
			g_error_free(error);
		}
	}
	else
	{
		/*
		 * search in [prefix]/share/[prgname] 
		 */
		char *bindir =
		   g_path_get_dirname(g_find_program_in_path(g_get_prgname()));
		char *filename = g_build_filename(bindir, "..", "share",
		   g_get_prgname(), gladeName,
		   NULL);
		if (g_file_test(filename, G_FILE_TEST_EXISTS))
		{
			main_window = gtk_builder_new();
			if (!gtk_builder_add_from_file(main_window, filename, &error))
			{
				g_warning("Couldn't load builder file: %s",
				   error->message);
				g_error_free(error);
			}
		}
		g_free(filename);
		g_free(bindir);
	}
	if (main_window == NULL)
		g_error("Unable to load Glade XML file for user interface");

	/*
	 * connect the signals in the interface 
	 */
	widget = GTK_WIDGET(gtk_builder_get_object(main_window, "button1"));
	g_signal_connect(widget, "clicked", G_CALLBACK(on_button1_clicked),
	   &p);

	widget = GTK_WIDGET(gtk_builder_get_object(main_window, "button2"));
	g_signal_connect(widget, "clicked", G_CALLBACK(on_button2_clicked),
	   NULL);

	GtkFileFilter *filter = gtk_file_filter_new();
	gtk_file_filter_add_pattern(filter, "*.txt");

	widget =
	   GTK_WIDGET(gtk_builder_get_object(main_window, "targetchooser"));
	gtk_file_chooser_add_filter((GtkFileChooser *) widget, filter);
	g_signal_connect(widget, "file-set",
	   G_CALLBACK(on_targetchooser_file_set), &p);

	widget =
	   GTK_WIDGET(gtk_builder_get_object(main_window, "arraychooser"));
	gtk_file_chooser_add_filter((GtkFileChooser *) widget, filter);
	g_signal_connect(widget, "file-set",
	   G_CALLBACK(on_arraychooser_file_set), &p);

	widget =
	   GTK_WIDGET(gtk_builder_get_object(main_window,
		  "instrumentchooser"));
	gtk_file_chooser_add_filter((GtkFileChooser *) widget, filter);
	g_signal_connect(widget, "file-set",
	   G_CALLBACK(on_instrumentchooser_file_set), &p);

	widget = GTK_WIDGET(gtk_builder_get_object(main_window, "SMchooser"));
	gtk_file_chooser_add_filter((GtkFileChooser *) widget, filter);
	g_signal_connect(widget, "file-set", G_CALLBACK(on_SMchooser_file_set),
	   &p);

	widget = GTK_WIDGET(gtk_builder_get_object(main_window, "HAchooser"));
	gtk_file_chooser_add_filter((GtkFileChooser *) widget, filter);
	g_signal_connect(widget, "file-set", G_CALLBACK(on_HAchooser_file_set),
	   &p);

	widget =
	   GTK_WIDGET(gtk_builder_get_object(main_window, "oifitschooser"));
	g_signal_connect(widget, "changed",
	   G_CALLBACK(on_oifitschooser_changed), &p);

	/*
	 * run the event loop 
	 */
	gtk_main();
}
