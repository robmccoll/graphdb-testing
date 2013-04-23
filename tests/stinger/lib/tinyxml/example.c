#include "tinyxml.h"
#include "tinyxml-dom.h"

#include <assert.h>
#include <stdio.h>

#define count(ARRAY) (sizeof (ARRAY) / sizeof *(ARRAY))

int main(void)
{
	struct tx_node nodes[32];
	char data[] = "<root><foo bar='spam'>test</foo><foo bar='eggs'/></root>";

	char *tail = tx_parse(data, count(nodes), nodes);
	assert(!tail); // check for parsing success

	{
	  puts("--- tx ---"); // use tx API
	  const char *path[] = { "root", "foo", "@bar" };
	  struct tx_node *current_bar = NULL;
	  while((current_bar = tx_get(nodes, current_bar, count(path), path)))
		  puts(current_bar->value);
	}
	{
	  puts("--- tx ---"); // use tx API
	  const char *path[] = { "root", "foo"};
	  struct tx_node *current_bar = NULL;
	  while((current_bar = tx_get(nodes, current_bar, count(path), path)))
		  puts(current_bar->value ? current_bar->value : "<NULL>");
	}

	puts("--- dom ---"); // use DOM wrapper
	dom_node_t foo_nodes[8];
	size_t count = dom_getElementsByTagName(
		nodes, "foo", count(foo_nodes), foo_nodes);

	size_t i = 0;
	for(; i < count; ++i)
		puts(dom_getAttribute(foo_nodes[i], "bar"));
}
